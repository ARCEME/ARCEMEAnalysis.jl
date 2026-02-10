module ARCEMEAnalysis
using Zarr: S3Store, Zarr
using Minio: MinioConfig
using YAXArrays: open_dataset, ⊘, xmap, XOutput, YAXArray, YAXArrays, Dataset
import DimensionalData as DD
using Colors: RGB
using Dates: DateTime, Year, Date
import CSV
using Statistics: mean
using DataStructures: SortedDict, counter
using ProgressMeter: @showprogress
using SpectralIndices: compute_index, SpectralIndices

include("download.jl")

const arceme_classes = SortedDict(
  0   => "No data",
  10  => "Tree cover",
  20  => "Shrubland",
  30  => "Grassland",
  40  => "Cropland",
  50  => "Built-up",
  60  => "Bare /sparse vegetation",
  70  => "Snow and Ice",
  80  => "Permanent water bodies",
  90  => "Herbaceous wetland",
  95  => "Mangroves",
  100 => "Moss and lichen",
)

export arceme_cubename, arceme_open, arceme_starttime, arceme_endtime, arceme_eventdate,
    arceme_coordinates, arceme_ndvi, arceme_rgb, arceme_eventlist, arceme_eventpairs, 
    arceme_classes, arceme_landcover, arceme_optical_band_fingerprints, arceme_radar_fingerprints,
    time_aggregate_fingerprint, arceme_validpairs, arceme_spectral

"""
    _arceme_cubenames(;batch="6")

List all available data cube names in the specified batch stored in the ARCEME S3 bucket.
"""
function _arceme_cubenames(; batch="ARCEME-DC-6")
    cubenames = if local_cubepath === nothing
        store = S3Store("$(batch)/", MinioConfig(httpstore))

        resp = Zarr.cloud_list_objects(store, batch)
        map(split(String(resp), "\n")) do p
            split(p, "/")[2]
        end
    else
        readdir(joinpath(local_cubepath, batch))
    end
    filter!(startswith("DC"), cubenames)
end


"""
    struct Event

Define a structure to hold information about an ARCEME event.
"""
struct Event
    uid::String
    dhp_label::Int
    dlabel::Int
    longitude::Float64
    latitude::Float64
    tp::Float64
    eventdate::DateTime
    source::Symbol
    unveg::Float64
    sparse::Float64
    teow::Int
    biome::String
end

"""
    arceme_eventlist()

Get a list of ARCEME events from the local CSV file.
"""
function arceme_eventlist()
    map(CSV.File(joinpath(@__DIR__, "..", "data", "dhp_global_subselection.csv"))) do row
        Event(
            row.uid,
            row.dhp_label,
            row.dlabel,
            row.longitude,
            row.latitude,
            row.tp_int,
            DateTime(row.date),
            Symbol(row.source),
            row.unveg,
            row.sparse,
            row.teow,
            row.biome,

        )
    end
end

arceme_cubename(event::Event) = 
"DC__$(event.uid)__$(Date(arceme_starttime(event)))__$(Date(arceme_endtime(event))).zarr"
# if event.source == :d
#     "DC__$(event.dlabel)_d__$(Date(arceme_starttime(event)))__$(Date(arceme_endtime(event))).zarr"
# else
#     "DC__$(event.dhp_label)_dhp__$(Date(arceme_starttime(event)))__$(Date(arceme_endtime(event))).zarr"
# end

"""
    arceme_eventpairs()

Get pairs of ARCEME events from the event list.
"""
arceme_eventpairs() =
    Iterators.partition(sort(arceme_eventlist(),by=i->(i.dhp_label,i.source)),2) |> collect

"""
    arceme_validpairs(;batch="ARCEME-DC-6")

Get valid pairs of ARCEME events from the local path (if set with arceme_set_localpath) or the http data store.
Default httpstore is "https://s3.waw3-2.cloudferro.com/swift/v1". It can be reset with arceme_set_httpstore.
"""
function arceme_validpairs(;batch="ARCEME-DC-6")
    allpairs = arceme_eventpairs()
    validpairs = if local_cubepath === nothing
        map(x -> all([Zarr.is_zgroup(Zarr.HTTPStore("$httpstore/$batch/$(arceme_cubename(i))"), "") for i in x]), allpairs)
    else
        map(x -> all([isfile(joinpath(local_cubepath, batch, string(arceme_cubename(i), ".zip"))) for i in x]), allpairs)
    end
    allpairs[validpairs]
end

function arceme_open(event::Event; batch="ARCEME-DC-6")
    arceme_open(arceme_cubename(event); batch=batch)
end

"""
     arceme_landcover(ds)

Get the land cover statistics for the ARCEME data cube dataset `ds`.
"""
function arceme_landcover(ds)
    cdr = counter(ds.ESA_LC);
    map(collect(arceme_classes)) do (k,v)
        count=get(cdr, k, 0)
        (key=k, class=v, count=count, fraction=count/1000000)
    end
end
arceme_landcover(ev::Event; batch="ARCEME-DC-6") = arceme_landcover(arceme_open(ev, batch=batch))


""" 
    arceme_open(cubename; batch="ARCEME-DC-6")

Open the specified ARCEME data cube from the local path (if with arceme_set_localpath) or 
    the httpstore (default "https://s3.waw3-2.cloudferro.com/swift/v1", reset with arceme_set_httpstore).
"""
function arceme_open(cubename; batch="ARCEME-DC-6")
    if local_cubepath === nothing
        open_dataset("$httpstore/$batch/$cubename", force_datetime=true)
    else
        open_dataset(joinpath(local_cubepath, batch, string(cubename, ".zip")))
    end
end


"""
    arceme_starttime(ds)

Get the start time of the ARCEME data cube dataset `ds`
"""
arceme_starttime(ds) = parse(DateTime, ds.properties["time_coverage_start"])
arceme_starttime(ev::Event) = ev.eventdate - Year(1)


"""
    arceme_endtime(ds)

Get the end time of the ARCEME data cube dataset `ds`
"""
arceme_endtime(ds) = parse(DateTime, ds.properties["time_coverage_end"])
arceme_endtime(ev::Event) = ev.eventdate + Year(1)

"""
    arceme_eventdate(ds)

Get the time of when the event in an ARCEME data cube dataset `ds` happened.
"""
arceme_eventdate(ds) = arceme_starttime(ds) + Year(1)
arceme_eventdate(ev::Event) = ev.eventdate

"""
    arceme_coordinates(ds)

Get the longitude and latitude coordinates of the center of the ARCEME data cube dataset `ds`.
"""
function arceme_coordinates(ds)
    trans = Proj.Transformation("EPSG:$(ds.properties["epsg"])", "OGC:84")
    x, y = ds.properties["central_x"], ds.properties["central_y"]
    lon, lat = round.(trans((x, y)), digits=4)
    return lon, lat
end


"""
    arceme_ndvi(ds)

Compute the NDVI (Normalized Difference Vegetation Index) for the ARCEME data cube dataset `ds`.
"""
function arceme_ndvi(ds)
    ndvi = broadcast(ds.B04, ds.B08, ds.cloud_mask, ds.SCL) do b4, b8, cl, scl
        (cl > 0 || (scl in (1, 3, 7, 8, 9, 10, 11))) && return NaN
        fb4 = boa(b4)
        fb8 = boa(b8)
        (fb8 - fb4) / (fb8 + fb4)
    end
    ds.cubes[:ndvi] = ndvi
    ds
end

#Helper functions to compute the indices from named tuples for type stability
compute_indexx(index::SpectralIndices.SpectralIndex{<:Any,B}, params::NamedTuple) where B = index.compute(Float64, params[B]...)
function listofindices(indices, values)
    map(indices) do index
        compute_indexx(index, values)
    end
end

"""
    arceme_spectral(ds, indices::Vector{String})

Compute the listed indices using SpectralIndices.jl. Not working.
"""
function arceme_spectral(ds, indices::Vector{String}; platform="sentinel2") 
    tmp = if platform == "sentinel2" || platform == "sentinel2a" || platform == "sentinel2b"
        pl = platform == "sentinel2" ? "sentinel2a" : platform
        _compute_indices(ds, indices, pl)
    elseif platform=="sentinel1"
        tmp = broadcast(ds.vv, ds.vh) do VV, VH
            compute_index(indices; VV,VH)
        end
    else
        error("platform $platform is not supported")
    end
    foreach(pairs(tmp)) do (k,v)
        ds.cubes[k] = v
    end
    ds
end

"""
    arceme_spectral(ds, index::String)

Compute index using SpectralIndices.jl, e.g. "NDVI". Not Working (when actually requesting the data).

Example:
validpairs = arceme_validpairs()
ds_d,ds_dhp = arceme_open.(validpairs[80])
arceme_spectral(ds_d, "NDVI")
@time ds_d.NDVI[x=1,y=1,].data[:]
ERROR: MethodError: no method matching (::XFunction{ARCEMEAnalysis.var"#36#37"{String}, XOutput{Tuple{}, Tuple{}, Int64}, Tuple{}})
The function `XFunction{ARCEMEAnalysis.var"#36#37"{String}, XOutput{Tuple{}, Tuple{}, Int64}, Tuple{}}(ARCEMEAnalysis.var"#36#37"{String}("NDVI"), XOutput{Tuple{}, Tuple{}, Int64}((), (), 1, Dict{Any, Any}()), (), false)` exists, but no method is defined for this combination of argument types.
"""
function arceme_spectral(ds, index::String; platform="sentinel2") 
    if platform=="sentinel2" || platform=="sentinel2a" || platform=="sentinel2b"
        tmp = broadcast(ds.cloud_mask, ds.SCL, ds.B01, ds.B02, ds.B03, ds.B04, ds.B08, ds.B05, ds.B06, ds.B07, ds.B11, ds.B12, ds.B09) do cl, scl, b1, b2, b3, b4, b8, b5, b6, b7, b11, b12, b9
            # BOA
            A = boa(b1); B = boa(b2); G = boa(b3); R = boa(b4); N = boa(b8)
            RE1 = boa(b5); RE2 = boa(b6); RE3 = boa(b7)
            S1 = boa(b11); S2 = boa(b12);  WV = boa(b9)
            # apply cloud masking here? (or after)
            (cl > 0 || (scl in (1, 3, 7, 8, 9, 10, 11))) && return NaN
            compute_index(index; A, B, G, R, N, RE1, RE2, RE3, S1, S2, WV, L=0.5)
        end
        ds.cubes[Symbol(index)] = tmp
    elseif platform=="sentinel1"
        tmp = broadcast(ds.vv, ds.vh) do VV, VH
            compute_index(index; VV,VH)
        end
    else
        error("platform $platform is not supported")
    end
    ds
end

"""
    boa(band; BOA_ADD_OFFSET = -1000, QUANTIFICATION_VALUE = 10000)

Compute the radiometric offsets for Sentinel 2 bands to get values at bottom of atmosphere.
The default values are valid for data processed from baseline 04.00 (January 2022) onwards.
The entire Sentinel-2 archive in CDSE has been reprocessed and is now available in baseline 05.xx, with consistent offset.
"""
boa(band; BOA_ADD_OFFSET = -1000, QUANTIFICATION_VALUE = 10000) = (band + BOA_ADD_OFFSET) / QUANTIFICATION_VALUE

"""
    arceme_rgb(ds)

Compute the RGB composite for the ARCEME data cube dataset `ds`.
"""
arceme_rgb(ds) =
    broadcast(ds.B02, ds.B03, ds.B04) do b, g, r
        # m = typemax(Int16)
        # RGB(r / m * 4, g / m * 4, b / m * 4)
        # RGB(min(1.0,r / m *4) , min(1.0,g / m *4) , min(1.0,b / m *4))
        RGB(boa(r), boa(g), boa(b))
end

"""
    arceme_optical_band_fingerprints(ds_d, ds_dhp)

Computes fingerprints for all optical bands in the ARCEME data cube datasets `ds_d` and `ds_dhp`. 
Also aggregates the fingerprints over time in 2-monthly windows for the whol 2 years around the event date.
"""
function arceme_optical_band_fingerprints(ds_d, ds_dhp)

    optical_bands = [:B01, :B02, :B03, :B04, :B05, :B06, :B07, :B08, :B8A, :B09, :B11, :B12]
    banddim = DD.Dim{:band}(string.(optical_bands))
    eventdate_d = arceme_eventdate(ds_d)
    eventdate_dhp = arceme_eventdate(ds_dhp)

    fp_d = @showprogress desc = "Drought fingerprint.." map(optical_bands) do band
        arceme_bias_corrected_fp(band, ds_d)
    end
    fingerprint_sparse_d = YAXArrays.concatenatecubes(map(i -> i.fp, fp_d), banddim)
    fingerprint_uncor_sparse_d = YAXArrays.concatenatecubes(map(i -> i.fp_uncorrected, fp_d), banddim)

    fp_dhp = @showprogress desc = "DHP fingerprint......" map(optical_bands) do band
        arceme_bias_corrected_fp(band, ds_dhp)
    end
    fingerprint_sparse_dhp = YAXArrays.concatenatecubes(map(i -> i.fp, fp_dhp), banddim)
    fingerprint_uncor_sparse_dhp = YAXArrays.concatenatecubes(map(i -> i.fp_uncorrected, fp_dhp), banddim)
    
    fingerprints_d = time_aggregate_fingerprint(fingerprint_sparse_d, eventdate_d, banddim, :time_sentinel_2_l2a)
    fingerprints_dhp = time_aggregate_fingerprint(fingerprint_sparse_dhp, eventdate_dhp, banddim, :time_sentinel_2_l2a)

    Dataset(; fingerprints_d, fingerprints_dhp, fingerprint_sparse_d, fingerprint_sparse_dhp, fingerprint_uncor_sparse_d, fingerprint_uncor_sparse_dhp)
end

"""
    arceme_radar_fingerprints(ds_d, ds_dhp)

Computes fingerprints for all radar bands in the ARCEME data cube datasets `ds_d` and `ds_dhp`.
Also aggregates the fingerprints over time in 2-monthly windows for the whol 2 years around the event date.
"""
function arceme_radar_fingerprints(ds_d, ds_dhp)

    radar_bands = [:vv, :vh]
    banddim = DD.Dim{:band}(string.(radar_bands))
    eventdate_d = arceme_eventdate(ds_d)
    eventdate_dhp = arceme_eventdate(ds_dhp)

    fp_d = map(radar_bands) do band
        arceme_uncorrected_fp(band, ds_d)
    end
    fingerprint_sparse_d = YAXArrays.concatenatecubes(map(i -> i.fp, fp_d), banddim)


    fp_dhp = map(radar_bands) do band
        arceme_uncorrected_fp(band, ds_dhp)
    end
    fingerprint_sparse_dhp = YAXArrays.concatenatecubes(map(i -> i.fp, fp_dhp), banddim)

    fingerprints_d = time_aggregate_fingerprint(fingerprint_sparse_d, eventdate_d, banddim, :time_sentinel_1_rtc)
    fingerprints_dhp = time_aggregate_fingerprint(fingerprint_sparse_dhp, eventdate_dhp, banddim, :time_sentinel_1_rtc)

    Dataset(; fingerprints_d, fingerprints_dhp, fingerprint_sparse_d, fingerprint_sparse_dhp)
end


"""
    arceme_merge_fingerprints(optical_fp, radar_fp)

Merges the optical and radar fingerprints datasets into a single dataset.
"""
function arceme_merge_fingerprints(optical_fp, radar_fp)
    newbands = DD.Dim{:band}([optical_fp.band.val; radar_fp.band.val])
    newarrays = map((:fingerprints_d,:fingerprints_dhp)) do cubename
        mergeddata = cat(optical_fp[cubename].data ./ typemax(UInt16), radar_fp[cubename].data, dims=3)
        cubename => YAXArray((radar_fp.lc, radar_fp.Ti, newbands), mergeddata)
    end
    Dataset(; newarrays...)
end

"""
    arceme_all_fingerprints(ds_d, ds_dhp)

Computes both optical and radar fingerprints for the ARCEME data cube datasets `ds_d` and `ds_dhp`.
Also aggregates the fingerprints over time in 2-monthly windows for the whol 2 years
around the event date.
"""
function arceme_all_fingerprints(ds_d, ds_dhp)
    optical_fp = arceme_optical_band_fingerprints(ds_d, ds_dhp)
    radar_fp = arceme_radar_fingerprints(ds_d, ds_dhp)
    arceme_merge_fingerprints(optical_fp, radar_fp)
end

function time_aggregate_fingerprint(allbands, eventdate, banddim, timeaxis)
    step_per_year = 6
    ts = Array(DD.dims(allbands, timeaxis))
    data = allbands.data
    timestep_to_group(t, eventdate, step_per_year) = (clamp(ceil(Int, (t - eventdate).value / (365.25 * 24 * 60 * 60 * 1000) * step_per_year), -step_per_year + 1, step_per_year))
    groupinds = timestep_to_group.(ts, eventdate, step_per_year)
    @show unique(groupinds)
    fingerprint_timesteps = (-step_per_year+1):(step_per_year)
    res = map(fingerprint_timesteps) do t
        ii = findall(==(t), groupinds)
        isnothing(ii) && return NaN
        mapslices(data[:, ii, :], dims=2) do sts
            all(ismissing, sts) && return NaN
            mean(skipmissing(sts))
        end
    end
    YAXArray((allbands.lc, DD.Ti((fingerprint_timesteps .- 0.5) ./ step_per_year .* 12), banddim), cat(res..., dims=2))
end



include("spatialdebias.jl")
include("spectral_helpers.jl")
end #module