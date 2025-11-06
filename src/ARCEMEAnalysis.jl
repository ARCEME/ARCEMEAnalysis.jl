module ARCEMEAnalysis
using Zarr: S3Store, Zarr
using Minio: MinioConfig
using YAXArrays: open_dataset, ⊘, xmap, XOutput, YAXArray
import YAXArrays
import DimensionalData as DD
using Colors: RGB
using Dates: DateTime, Year, Date
import CSV
using DataStructures: SortedDict, counter
import Proj
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
    arceme_classes, arceme_landcover

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
    arceme_open(cubename; batch="SECONDBATCH")

Open the specified ARCEME data cube from the S3 bucket.
"""
arceme_open(cubename; batch="SECONDBATCH") =
    open_dataset("https://s3.waw3-2.cloudferro.com/swift/v1/ARCEME-DATACUBES/$batch/$cubename")


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
arceme_ndvi(ds) =
    broadcast(ds.B04, ds.B08, ds.cloud_mask) do b4, b8, cl
        cl > 0 && return NaN
        fb4 = b4 / typemax(Int16)
        fb8 = b8 / typemax(Int16)
        (fb8 - fb4) / (fb8 + fb4)
    end

"""
    arceme_rgb(ds)

Compute the RGB composite for the ARCEME data cube dataset `ds`.
"""
arceme_rgb(ds) =
    broadcast(ds.B02, ds.B03, ds.B04) do b, g, r
        m = typemax(Int16)
        RGB(r / m * 8, g / m * 8, b / m * 8)
end

include("spatialdebias.jl")
end #module