module ARCEMEAnalysis
using Zarr: S3Store, Zarr
using Minio: MinioConfig
using YAXArrays: open_dataset
using Colors: RGB
using Dates: DateTime, Year

export arceme_cubenames, arceme_open, arceme_starttime, arceme_endtime, arceme_eventdate, arceme_coordinates, arceme_ndvi, arceme_rgb

"""
    arceme_cubenames(;batch="SECONDBATCH")

List all available data cube names in the specified batch stored in the ARCEME S3 bucket.
"""
function arceme_cubenames(; batch="SECONDBATCH")
    store = S3Store("ARCEME-DATACUBES/", MinioConfig("https://s3.waw3-2.cloudferro.com/swift/v1"))

    resp = Zarr.cloud_list_objects(store, batch)
    cubenames = map(split(String(resp), "\n")) do p
        split(p, "/")[2]
    end
    filter!(startswith("DC"), cubenames)
end

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

"""
    arceme_endtime(ds)

Get the end time of the ARCEME data cube dataset `ds`
"""
arceme_endtime(ds) = parse(DateTime, ds.properties["time_coverage_end"])

"""
    arceme_eventdate(ds)

Get the time of when the event in an ARCEME data cube dataset `ds` happened.
"""
arceme_eventdate(ds) = arceme_starttime(ds) + Year(1)


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
        RGB(r / m * 4, g / m * 4, b / m * 4)
end


end #module