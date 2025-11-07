using ARCEMEAnalysis
using Test

@testset "ARCEMEAnalysis.jl" begin
    # Write your tests here.
    eventlist = arceme_eventlist()
    eventpairs = arceme_eventpairs()
    cubename = arceme_cubename(eventpairs[105][2])
    ds = arceme_open(cubename)
    ndvi = arceme_ndvi(ds)
    lccube = arceme_landcover(ds)
    cloudcube = ds.cloud_mask
    @time fp = arceme_bias_corrected_fp(ndvi, cloudcube,lccube)
end
