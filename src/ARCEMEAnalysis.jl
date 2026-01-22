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
    arceme_classes, arceme_landcover, arceme_optical_band_footprints, arceme_radar_footprints,
    time_aggregate_footprint

"""
    _arceme_cubenames(;batch="SECONDBATCH")

List all available data cube names in the specified batch stored in the ARCEME S3 bucket.
"""
function _arceme_cubenames(; batch="SECONDBATCH")
    store = S3Store("ARCEME-DATACUBES/", MinioConfig("https://s3.waw3-2.cloudferro.com/swift/v1"))

    resp = Zarr.cloud_list_objects(store, batch)
    cubenames = map(split(String(resp), "\n")) do p
        split(p, "/")[2]
    end
    filter!(startswith("DC"), cubenames)
end


"""
    struct Event

Define a structure to hold information about an ARCEME event.
"""
struct Event
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
if event.source == :d
    "DC__$(event.dlabel)_d__$(Date(arceme_starttime(event)))__$(Date(arceme_endtime(event))).zarr"
else
    "DC__$(event.dhp_label)_dhp__$(Date(arceme_starttime(event)))__$(Date(arceme_endtime(event))).zarr"
end

"""
    arceme_eventpairs()

Get pairs of ARCEME events from the event list.
"""
arceme_eventpairs() =
    Iterators.partition(sort(arceme_eventlist(),by=i->(i.dhp_label,i.source)),2) |> collect

function arceme_open(event::Event)
    arceme_open(arceme_cubename(event))
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
arceme_landcover(ev::Event) = arceme_landcover(arceme_open(ev))


""" 
    arceme_open(cubename; batch="ARCEME-DATACUBES/SECONDBATCH")

Open the specified ARCEME data cube from the S3 bucket.
"""
arceme_open(cubename; batch="ARCEME-DC-5") =
    open_dataset("https://s3.waw3-2.cloudferro.com/swift/v1/$batch/$cubename", force_datetime=true)


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
        fb4 = b4 / typemax(Int16)
        fb8 = b8 / typemax(Int16)
        (fb8 - fb4) / (fb8 + fb4)
    end
    oldcubes = ds.cubes
    ds.cubes[:ndvi] = ndvi
    ds
end

"""
    arceme_rgb(ds)

Compute the RGB composite for the ARCEME data cube dataset `ds`.
"""
arceme_rgb(ds) =
    broadcast(ds.B02, ds.B03, ds.B04) do b, g, r
        m = typemax(Int16)
        # RGB(r / m * 4, g / m * 4, b / m * 4)
        # RGB(min(1.0,r / m *4) , min(1.0,g / m *4) , min(1.0,b / m *4))
        RGB(r / m , g / m , b / m )
end

"""
    arceme_optical_band_footprints(ds_d, ds_dhp)

Computes footprints for all optical bands in the ARCEME data cube datasets `ds_d` and `ds_dhp`. 
Also aggregates the footprints over time in 2-monthly windows for the whol 2 years around the event date.
"""
function arceme_optical_band_footprints(ds_d, ds_dhp)

    optical_bands = [:B01, :B02, :B03, :B04, :B05, :B06, :B07, :B08, :B8A, :B09, :B11, :B12]
    banddim = DD.Dim{:band}(string.(optical_bands))
    eventdate_d = arceme_eventdate(ds_d)
    eventdate_dhp = arceme_eventdate(ds_dhp)

    fp_d = @showprogress desc = "Drought Footprint.." map(optical_bands) do band
        arceme_bias_corrected_fp(band, ds_d)
    end
    footprint_sparse_d = YAXArrays.concatenatecubes(map(i -> i.fp, fp_d), banddim)
    footprint_uncor_sparse_d = YAXArrays.concatenatecubes(map(i -> i.fp_uncorrected, fp_d), banddim)

    fp_dhp = @showprogress desc = "DHP Footprint......" map(optical_bands) do band
        arceme_bias_corrected_fp(band, ds_dhp)
    end
    footprint_sparse_dhp = YAXArrays.concatenatecubes(map(i -> i.fp, fp_dhp), banddim)
    footprint_uncor_sparse_dhp = YAXArrays.concatenatecubes(map(i -> i.fp_uncorrected, fp_dhp), banddim)
    
    footprints_d = time_aggregate_footprint(footprint_sparse_d, eventdate_d, banddim, :time_sentinel_2_l2a)
    footprints_dhp = time_aggregate_footprint(footprint_sparse_dhp, eventdate_dhp, banddim, :time_sentinel_2_l2a)

    Dataset(; footprints_d, footprints_dhp, footprint_sparse_d, footprint_sparse_dhp, footprint_uncor_sparse_d, footprint_uncor_sparse_dhp)
end

"""
    arceme_radar_footprints(ds_d, ds_dhp)

Computes footprints for all radar bands in the ARCEME data cube datasets `ds_d` and `ds_dhp`.
Also aggregates the footprints over time in 2-monthly windows for the whol 2 years around the event date.
"""
function arceme_radar_footprints(ds_d, ds_dhp)

    radar_bands = [:vv, :vh]
    banddim = DD.Dim{:band}(string.(radar_bands))
    eventdate_d = arceme_eventdate(ds_d)
    eventdate_dhp = arceme_eventdate(ds_dhp)

    fp_d = map(radar_bands) do band
        arceme_uncorrected_fp(band, ds_d)
    end
    footprint_sparse_d = YAXArrays.concatenatecubes(map(i -> i.fp, fp_d), banddim)


    fp_dhp = map(radar_bands) do band
        arceme_uncorrected_fp(band, ds_dhp)
    end
    footprint_sparse_dhp = YAXArrays.concatenatecubes(map(i -> i.fp, fp_dhp), banddim)

    footprints_d = time_aggregate_footprint(footprint_sparse_d, eventdate_d, banddim, :time_sentinel_1_rtc)
    footprints_dhp = time_aggregate_footprint(footprint_sparse_dhp, eventdate_dhp, banddim, :time_sentinel_1_rtc)

    Dataset(; footprints_d, footprints_dhp, footprint_sparse_d, footprint_sparse_dhp)
end


"""
    arceme_merge_footprints(optical_fp, radar_fp)

Merges the optical and radar footprints datasets into a single dataset.
"""
function arceme_merge_footprints(optical_fp, radar_fp)
    newbands = DD.Dim{:band}([optical_fp.band.val; radar_fp.band.val])
    newarrays = map((:footprints_d,:footprints_dhp)) do cubename
        mergeddata = cat(optical_fp[cubename].data ./ typemax(UInt16), radar_fp[cubename].data, dims=3)
        cubename => YAXArray((radar_fp.lc, radar_fp.Ti, newbands), mergeddata)
    end
    Dataset(; newarrays...)
end

"""
    arceme_all_footprints(ds_d, ds_dhp)

Computes both optical and radar footprints for the ARCEME data cube datasets `ds_d` and `ds_dhp`.
Also aggregates the footprints over time in 2-monthly windows for the whol 2 years
around the event date.
"""
function arceme_all_footprints(ds_d, ds_dhp)
    optical_fp = arceme_optical_band_footprints(ds_d, ds_dhp)
    radar_fp = arceme_radar_footprints(ds_d, ds_dhp)
    arceme_merge_footprints(optical_fp, radar_fp)
end

function time_aggregate_footprint(allbands, eventdate, banddim, timeaxis)
    step_per_year = 6
    ts = Array(DD.dims(allbands, timeaxis))
    data = allbands.data
    timestep_to_group(t, eventdate, step_per_year) = (clamp(ceil(Int, (t - eventdate).value / (365.25 * 24 * 60 * 60 * 1000) * step_per_year), -step_per_year + 1, step_per_year))
    groupinds = timestep_to_group.(ts, eventdate, step_per_year)
    @show unique(groupinds)
    footprint_timesteps = (-step_per_year+1):(step_per_year)
    res = map(footprint_timesteps) do t
        ii = findall(==(t), groupinds)
        isnothing(ii) && return NaN
        mapslices(data[:, ii, :], dims=2) do sts
            all(ismissing, sts) && return NaN
            mean(skipmissing(sts))
        end
    end
    YAXArray((allbands.lc, DD.Ti((footprint_timesteps .- 0.5) ./ step_per_year .* 12), banddim), cat(res..., dims=2))
end



include("spatialdebias.jl")
end #module