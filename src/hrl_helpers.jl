module HRL
using ARCEMEAnalysis
import CondaPkg
using  PythonCall: pyimport, pyconvert, pylist
import GeoJSON
using Preferences: @set_preferences!, @load_preference
import ArchGDAL as AG
using YAXArrays
using Dates
using ProgressMeter
using DataStructures: SortedDict, counter

"""
hrl_set_localpath(path)

Sets the local path to download HRL data.
If no local preference is set, the default location is used: pwd()

Setting this preference needs recompilation.
"""
function hrl_set_localpath(path)
    @set_preferences!("hrl_localpath" => path)
    @warn "Preferences changed! Restart Julia for this change to take effect."
end

const hrl_localpath = @load_preference("hrl_localpath", pwd())

if !isdir(hrl_localpath)
  mkdir(hrl_localpath)
end


"""
`hrl_hda()`

Sets up client to Harmonised Data Access API on WEkEO
Returns HDA client as Python Object
"""
function hrl_hda()
# Harmonised Data Access API on WEkEO
  hda = try pyimport("hda") #("Client", "Configuration")
  catch
      CondaPkg.add("hda")
      pyimport("hda")
  end
  # Configure user's credentials with a .hdarc file
  hda_client = try hda.Client()
  catch
      @error "Please configure your WEkEO user's credentials with a .hdarc file  (see https://hda.readthedocs.io/en/latest/usage.html#client-configuration for details)"
  end
  return hda_client
end

"""
`hda_datasets(;filteron="HRL")`

Search for all datasets available on WEkEO in whose name `filteron` occurs. Very slow!

Returns a `Vector` of `String`s
"""
function hda_datasets(;filteron="HRL")
  hda_client = hrl_hda()
  datasets = hda_client.datasets() 
  dataset_ids = [item["dataset_id"] for item in datasets] 
  # convert to julia
  pyconvert(Vector{String}, dataset_ids) |> filter(x -> occursin(filteron, x))
end

"""
`hda_products(dataset_id)`

Lists all valid productType for `dataset_id` on WEkEO.
Returns a vector of strings
"""
function hda_products(dataset_id)
  hda_client = hrl_hda()
  metadata = hda_client.metadata(dataset_id)
  map(x -> x["const"], pyconvert(Vector{Dict},metadata["properties"]["productType"]["oneOf"]))
end

"""
`hda_years(dataset_id)`

Lists all valid years or epochs for `dataset_id` on WEkEO.
Returns a vector of strings
"""
function  hda_years(dataset_id)
  hda_client = hrl_hda()
  metadata = hda_client.metadata(dataset_id)
  map(x -> x["const"], pyconvert(Vector{Dict},metadata["properties"]["year"]["oneOf"]))
end


"""
`hrl_download(dataset_id::String, coordinates::Tuple, yr::Union{String,Int}; extent=0.05, product::{Union{Nothing,String}}=nothing)`
`hrl_download(dataset_id::Vector{String}, coordinates::Tuple, yr::Union{String,Int}; extent=0.05, product::{Union{Nothing,String}}=nothing)`

Downloads Copernicus Land Monitoring Service High Resolution Layer products from WEkEO using the Harmonised Data Access API. 
  `dataset_id` must be a single valid dataset (e.g., Cropland: "EO:EEA:DAT:HRL:CRL"; Tree cover and Forest: "EO:EEA:DAT:HRL:TCF"; 
  Grassland: "EO:EEA:DAT:HRL:GRA") or a vector of vallid datasets.
  Search for all tiles near `coordinates` with (lon, lat) +- `extent`, for year `yr`. Some products are only available for three-years epochs.
  Data are downloaded at `hrl_localpath` only if they do not already exist.
  By default, all products from `dataset_id` are returned. To search for a specific product, set `product` to a valid productType (e.g. "Crop Types")
  
Returns a Vector of filenames.

See also `hda_datasets`, `hda_products`, `hda_years`
"""
function hrl_download(dataset_id::String, coord::Tuple, yr::Union{String,Int}; extent=0.05, product::Union{Nothing,String}=nothing)
    hda_client = hrl_hda()
    json = pyimport("json")

    (lon,lat) = coord
    fnames = []

    # Build query as a JSON object
    if !isnothing(product)
      pt = """,
      "productType": "$(product)"
      """
    else
      pt = ""
    end

    query = """
    {
      "dataset_id": "$(dataset_id)",
      "bbox": [
        $(lon-extent),
        $(lat-extent),
        $(lon+extent),
        $(lat+extent)
      ],
      "resolution": "10m",
      "year": "$yr",
      "itemsPerPage": 200,
      "startIndex": 0$(pt)
    }
    """
    
    matches = hda_client.search(json.loads(query))
    try
    for item in matches.results
      # println(item["id"])
      filename = pyconvert(String, item["id"])*".zip"
      # store info:
      push!(fnames, filename)
    end
    catch
      @warn "No matches found for $(dataset_id) at $(coordinates) in year $(yr)"
      return fnames
    end

    # Download
    todownload = [!isfile(joinpath(hrl_localpath, pyconvert(String, item["id"])*".zip")) for  item in matches.results]
    @showprogress desc="Downloading..." for i in pylist(findall(todownload).-1)
    #   println("Downloading $(matches.results[i]["id"])")
      matches[i].download(download_dir=hrl_localpath) 
    end
    return fnames
end

function hrl_download(datasets::Vector{String}, coord::Tuple, yr::Union{String,Int}; kwargs...)
  fnames=[];
  for dataset_id in datasets
      append!(fnames, HRL.hrl_download(dataset_id, coord, yr; kwargs...))
  end
  return fnames
end

const hrl_dict = Dict(
    "CTY" => "Crop types",
    "CTYCL" => "Crop types confidence level",
    "CPMCH" => "Main Crop Harvest",
    "CPMCHCL" => "Main Crop Harvest confidence level",
    "CPMCE" => "Main Crop Emergence",
    "CPMCECL" => "Main Crop Emergence confidence level",
    "CPMCD" => "Main Crop Duration",
    "CPMCDCL" => "Main Crop Duration confidence level",
    "CPBSB" => "Bare soil before",
    "CPBSBCL" => "Bare soil before confidence level",
    "CPBSA" => "Bare soil after",
    "CPBSACL" => "Bare soil after confidence level",
    "CPSCT" => "Secondary Crop Type",
    "CPSCE" => "Secondary Crop Emergence",
    "CPSCD" => "Secondary Crop Duration",
    "CPSCDCL" => "Secondary Crop Duration confidence level",
    "CPFLP" => "Fallow Land Presence",
    "CPFLPCL" => "Fallow Land Presence confidence level",
    # "CPFLD" => "Fallow Land Duration",
    # "CPFLDCL" => "Fallow Land Duration confidence level",
    "CPCSY" => "Cropping Seasons Yearly",
    # Forest
    "DLT" => "Dominant Leaf Type",
    "DLTCL" => "Dominant Leaf Type Confidence Layer",
    "TCD" => "Tree Cover Density",
    "TCDCL" => "Tree Cover Density Confidence Layer",
    # Grassland
    "GRA" => "Grassland",
    "GRACL" => "Grassland Confidence Layer",
    "HER" => "Herbaceous Cover",
    "PLOUGH" => "Ploughing Indicator",
    "GRAME"  => "Grassland Mowing Events",
    "GRAMECL" => "Grassland Mowing Events Confidence Layer",
    "GRAMD1" => "Grassland Mowing Dates of first mowing event",
    "GRAMD2" => "Grassland Mowing Dates of second mowing event",
    "GRAMD3" => "Grassland Mowing Dates of third mowing event",
    "GRAMD4" => "Grassland Mowing Dates of fourth mowing event",
)

const hrl_legends = Dict(
    "MCTY" => SortedDict(
        1 => "No cropland",
        2 => "Annual crop",
        3 => "Permanent crop"
    ),
    "CTY" => SortedDict(
        0 =>  "No cropland",
        1110 => "Wheat",
        1120 => "Barley",
        1130 => "Maize",
        1140 => "Rice",
        1150 => "Other cereals",
        1210 => "Fresh Vegetables",
        1220 => "Dry Pulses",
        1310 => "Potatoes",
        1320 => "Sugar Beet",
        1410 => "Sunflower",
        1420 => "Soybeans",
        1430 => "Rapeseed",
        1440 => "Flax cotton and hemp",
        2100 => "Grapes",
        2200 => "Olives",
        2310 => "Fruits",
        2320 => "Nuts",
        3100 => "Unclassified annual crop",
        3200 => "Unclassified permanent crop",
        65535 => "outside area"
    ),
    "CTYCL" => "0-100 Probability expressed as a percentage. 253 No Cropland. 255 Outside area",
    "CPMCH" => "Harvest date as YYDOY where YY = last 2 digits of the year (e.g., 19 for 2019) and DOY is the day of the year (1-366). 
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        65535 Outside area.",
    "CPMCHCL" => "Uncertainty in days (1-40). 
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        64534 No confidence could be calculated.
        65535 Outside area.",   
    "CPMCE" => "Emergence date as YYDOY where YY = last 2 digits of the year (e.g., 19 for 2019) and DOY is the day of the year (1-366). 
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        65535 Outside area.",
    "CPMCECL" => "Uncertainty in days (1-40). 
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        64534 No confidence could be calculated.
        65535 Outside area.",
    "CPMCD" => "Days (40-366).  
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        65535 Outside area.",
    "CPMCDCL" => "Uncertainty in days (1-80). 
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        64534 No confidence could be calculated.
        65535 Outside area.",
    "CPBSB" => "Days (1-295).  
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65529 Period outside calendar boundaries. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        65535 Outside area.",
    "CPBSBCL" => "Uncertainty in days (1-80). 
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65529 Period outside calendar boundaries. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        64534 No confidence could be calculated.
        65535 Outside area.",
    "CPBSA" => "Days (1-295).  
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65529 Period outside calendar boundaries. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        65535 Outside area.",
    "CPBSACL" => "Uncertainty in days (1-80). 
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65529 Period outside calendar boundaries. 
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        64534 No confidence could be calculated.
        65535 Outside area.",
    "CPSCT" => SortedDict(
        0 => "No annual cropland", 
        1 => "Short Summer",
        2 => "Long Summer",
        3 => "Short Winter",
        4 => "Long Winter", 
        65526 => "Fallow land", 
        65527 => "No cropping pattern detected",  
        65530 => "No secondary crop growing season delineated",
        65531 => "Not enough data", 
        65532 => "No cropping season detected",  
        65533 => "Growing season extends beyond timeframe",
        65535 => "Outside area"
    ),
    "CPSCE" => "Emergence date as YYDOY where YY = last 2 digits of the year (e.g., 19 for 2019) and DOY is the day of the year (1-366). 
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65530 No secondary crop growing season delineated.
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        65535 Outside area.",
    "CPSCD" => "Days (40-366).  
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65530 No secondary crop growing season delineated.
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        65535 Outside area.",
    "CPSCDCL" => "Uncertainty in days (1-80). 
        0 No annual cropland. 
        65526 Fallow land. 
        65527 No cropping pattern detected. 
        65530 No secondary crop growing season delineated.
        65531 Not enough data. 
        65532 No cropping season detected. 
        65533 Growing season extends beyond timeframe. 
        64534 No confidence could be calculated.
        65535 Outside area.",
    "CPFLP" => SortedDict(
        0 => "No fallow land",
        1 => "Fallow land",
        65535 => "Outside area"
    ),
    "CPFLPCL" => "0-100 Probability expressed as a percentage. 253 No Cropland. 65535 Outside area",
    "CPCSY" => SortedDict(
        0 => "No annual cropland", 
        1 => "One growing season",
        2 => "Two growing seasons",
        65526 => "Fallow land", 
        65527 => "No cropping pattern detected",  
        65531 => "Not enough data", 
        65532 => "No cropping season detected",  
        65533 => "Growing season extends beyond timeframe",
        65535 => "Outside area"
    ),
    # Forest
    "DLT" => SortedDict(
        0 => "all non-tree covered areas",
        1 => "broadleaved trees",
        2 => "coniferous trees",
        255 => "outside area",
    ),
    "DLTCL" => "0-100: classification confidence. 
            253: all non-tree covered areas. 
            255: outside area",
    "TCD" => "0: all non-tree covered areas. 
            1-100: tree cover density in %. 
            255: outside area",
    "TCDCL" => "0-100: standard deviation of TCD estimate. 
            253: all non-tree covered areas
            255: outside area",
    # Grassland
    "GRA" => SortedDict(
        0 => "all non-grassland areas",
        1 => "grassland", 
        255 => "outside area",
    ),
    "GRACL" => "0-100: Classification confidence. 253: All non-grassland areas. 255: outside area",
    "HER" => SortedDict(
        0 => "non-grassland in reference year",
        1 => "temporary grassland in reference year",
        255 => "outside area",
    ),
    "PLOUGH" => "0: Indication of ploughing in current year. 
    1-6: Number of years since last indication of ploughing. 
    100: Change in herbaceous cover. 
    253: no ploughing information. 
    255: outside area",
    "GRAME"  => SortedDict(
        0 => "no mowing detected", 
        1 => "1 mowing event detected",
        2 => "2 mowing events",
        3 => "3 mowing events",
        4 => "4 mowing events",
        253 => "all non-herbaceous areas",
        255 => "outside area"
    ),
    "GRAMECL" => "0-100: Mowing detection confidence. 253: All non-mowing areas. 255: outside area",
    "GRAMD1" => "0: no mowing detected. 1-366: Start (DOY) of each mowing event. 65533: all non-herbaceous areas. 65535: outside area",
    "GRAMD2" => "0: no mowing detected. 1-366: Start (DOY) of each mowing event. 65533: all non-herbaceous areas. 65535: outside area",
    "GRAMD3" => "0: no mowing detected. 1-366: Start (DOY) of each mowing event. 65533: all non-herbaceous areas. 65535: outside area",
    "GRAMD4" => "0: no mowing detected. 1-366: Start (DOY) of each mowing event. 65533: all non-herbaceous areas. 65535: outside area",
)

"""
`hrl_warp(event::ARCEMEAnalysis.Event; batch="ARCEME-DC-6", dataset_id::Union{String, Vector{String}}=["EO:EEA:DAT:HRL:CRL", "EO:EEA:DAT:HRL:TCF", "EO:EEA:DAT:HRL:GRA"])`
`hrl_warp(cubename, fnames::Vector{String}; batch="ARCEME-DC-6")`
`hrl_warp(cubename::String; batch="ARCEME-DC-6", dataset_id::Union{String, Vector{String}}=["EO:EEA:DAT:HRL:CRL", "EO:EEA:DAT:HRL:TCF", "EO:EEA:DAT:HRL:GRA"])`

reprojects and clips files in `fnames` to match `event`'s coordinate reference system and extent,
 and saves them as a Zarr group in `batch-HRL` under `ARCEMEAnalysis.local_cubepath`.
"""
function hrl_warp(cubename, fnames::Vector; batch="ARCEME-DC-6")
  ds = arceme_open(cubename; batch)
  event_date = event_date = arceme_eventdate(ds)
  tmp = tempname(;suffix=".tif")
  res = ds.properties["resolution"]
  tmp = tempname(;suffix=".tif")

  for hrl in keys(hrl_dict)
    hrl_fnames = fnames[map(contains("$(hrl)_"), fnames)]
    if isempty(hrl_fnames)
        continue
    end
    run(Cmd(
        `gdalwarp 
            -t_srs $(AG.toPROJ4(AG.importEPSG(ds.properties["epsg"]))) 
            -tr $res $res
            -te $(minimum(ds.x)) $(minimum(ds.y)) $(maximum(ds.x)+res) $(maximum(ds.y)+res) 
            -r near
            $(["/vsizip/$(joinpath(hrl_localpath, hrl_fnames[i], "$(split(hrl_fnames[i],".")[1]).tif") )" for i in eachindex(hrl_fnames)]) 
            $tmp 
            -overwrite 
            `))
    c = Cube(tmp)
    # write to zarr
    props = Dict(
        "long_name" => hrl_dict["$(split(hrl_fnames[1],"_")[3])"],
        "year_coverage" => year(event_date),
        "legend" => hrl_legends["$(split(hrl_fnames[1],"_")[3])"],
        "epsg" => ds.properties["epsg"], 
        "resolution" => 10,
        "unit" => "meter",
        "source" => "$(join(hrl_fnames, ", ")) © European Union, Copernicus Land Monitoring Service 2025, European Environment Agency (EEA)",
        "processing" => "Downloaded from wEkEO. Reprojected to EPSG $(ds.properties["epsg"]) using nearest neighour, mosaicked and clipped to AOI with gdalwarp",
        "processed_by" => "M. Weynants, MPI-BGC",
        "processed_date" => string(Date(now())),
    )
    c2 = YAXArray((ds.x,ds.y), c.data, props, )
    nt = (Symbol("$(split(hrl_fnames[1],"_")[3])")=>c2,)
    tds = Dataset(;nt...)
    tds = setchunks(tds, (500,500))
    savedataset(tds, 
        path=joinpath(ARCEMEAnalysis.local_cubepath, "$(batch)-HRL", "$(split(cubename,".")[1]).zarr"), 
        append=true)
    end
    # zip
    run(Cmd(`zip -0 -r ../$(split(cubename,".")[1]).zarr.zip .`, dir=joinpath(ARCEMEAnalysis.local_cubepath, "$(batch)-HRL", "$(split(cubename,".")[1]).zarr")))
    rm(joinpath(ARCEMEAnalysis.local_cubepath, "$(batch)-HRL", "$(split(cubename,".")[1]).zarr"), recursive=true)
    return joinpath(ARCEMEAnalysis.local_cubepath, "$(batch)-HRL", "$(split(cubename,".")[1]).zarr.zip")
end

hrl_warp(event::ARCEMEAnalysis.Event; kwargs...) = hrl_warp(arceme_cubename(event); kwargs...)

function hrl_warp(cubename::String; batch="ARCEME-DC-6", dataset_id::Union{String, Vector{String}}=["EO:EEA:DAT:HRL:CRL", "EO:EEA:DAT:HRL:TCF", "EO:EEA:DAT:HRL:GRA"])
  ds = arceme_open(cubename; batch=batch)
  coord = arceme_coordinates(ds)
  yr = year(arceme_eventdate(ds))
  fnames = hrl_download(dataset_id, coord, yr)
  hrl_warp(cubename, fnames; batch)
end

"""
     arceme_stats(ds::YAXArrays.Dataset, name; classes=nothing)
     arceme_stats(ev::Event, name; batch="ARCEME-DC-6", classes=nothing)

Get the  statistics of layer `name` for the ARCEME data cube `ds`.
Default classes are extracted from `hrl_legends[name]`. Alternative classes can be provided as a dictionary.
"""
function arceme_stats(ds, name; classes=nothing)
    if isnothing(classes)
        classes = hrl_legends[name]
    end
    cdr = counter(ds[name]);
    map(collect(classes)) do (k,v)
        count=get(cdr, k, 0)
        (key=k, class=v, count=count, fraction=count/1000000)
    end
end
arceme_stats(ev::ARCEMEAnalysis.Event, name; batch="ARCEME-DC-6", classes=nothing) = arceme_stats(arceme_open(ev, batch=batch), name; classes = classes)

function ctykeymap(k)
    k == 0 && return 1
    k > 3000 && return 1
    ctL1 = k ÷ 1000
    ctL2 = (k - 1000*ctL1) ÷ 100
    ctL3 = (k - 1000*ctL1 - 100*ctL2) ÷ 10
    ctL1==1 && ctL2==1 && return ctL3+1 # cereals
    ctL1==1 && ctL2==2 && return ctL3+6 # vegetables and pulses
    ctL1==1 && ctL2==3 && return ctL3+8 # roots (potato, sugar beet)
    ctL1==1 && ctL2==4 && return ctL3+10 # oil (sunflower, soy, rapeseed)
    ctL1==2 && ctL2==1 && return ctL3+15 # fiber (flax, cotton, hemp)
    ctL1==2 && ctL2==2 && return ctL2+14 # grape, olive
    ctL1==2 && ctL2==3 && return ctL3+16 # fruit, nuts
end

function mctykeymap(k)
    k == 0 && return 1 # non crop
    k == 65535 && return 1
    ctL1 = k ÷ 1000
    ctL2 = (k - 1000*ctL1) ÷ 100
    ctL3 = (k - 1000*ctL1 - 100*ctL2) ÷ 10
    ctL1 == 1 && return 2 # annual crops
    # ctL1 == 2 && ctL2 == 1 && return 2 # annual fiber
    ctL1 == 2 && return 3 # permanent crop
    k == 3100 && return 2 # Unclassified annual crop
    k == 3200 && return 3 # Unclassified permanent crop
end
export hrl_download, hrl_set_localpath, hda_datasets, hda_products, hda_years, hrl_warp, arceme_stats
end #module