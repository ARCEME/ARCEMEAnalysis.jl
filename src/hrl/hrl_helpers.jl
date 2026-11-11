module HRL
using ARCEMEAnalysis
import CondaPkg
using PythonCall: pyimport, pyconvert, pylist
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
Mapping from old WEkEO HRL dataset/product identifiers to new CLMS dataset IDs.

Key values are `Vector{String}` — either a single-element vector for a direct product
mapping, or a multi-element vector for old parent datasets (CRL, TCF, GRA) that
returned all their products.

Old parent datasets resolved to vectors of all their products:
- `EO:EEA:DAT:HRL:CRL` → all cropland CLMS datasets (22 products, see hrl_dict lines
  before `#Forest`)
- `EO:EEA:DAT:HRL:TCF` → tree cover & forest CLMS datasets (TCD, TCDCL, DLT, DLTCL)
- `EO:EEA:DAT:HRL:GRA` → grassland CLMS datasets (GRA, GRACL, HER, PLOUGH, GRAME,
  GRAMECL, GRAMD1-4)

Confidence layers without corresponding CLMS datasets map to `String[]` (empty).
"""
const hrl_dataset_map = Dict{String,Vector{String}}(
    # --- Cropland parent (EO:EEA:DAT:HRL:CRL) ---
    # Returns ALL cropland products for backward compatibility
    "EO:EEA:DAT:HRL:CRL" => [
        "EO:CLMS:DAT:CTY_10M_YEARLY_2017_PRESENT",          # CTY
        "",                                                   # CTYCL → no CLMS dataset (confidence)
        "EO:CLMS:DAT:CPMCH_10M_YEARLY_2017_PRESENT",        # CPMCH
        "",                                                   # CPMCHCL → no CLMS dataset
        "EO:CLMS:DAT:CPMCE_10M_YEARLY_2017_PRESENT",        # CPMCE
        "",                                                   # CPMCECL → no CLMS dataset
        "EO:CLMS:DAT:CPMCD_10M_YEARLY_2017_PRESENT",        # CPMCD
        "",                                                   # CPMCDCL → no CLMS dataset
        "EO:CLMS:DAT:CPBSB_10M_YEARLY_2017_PRESENT",        # CPBSB
        "",                                                   # CPBSBCL → no CLMS dataset
        "EO:CLMS:DAT:CPBSA_10M_YEARLY_2017_PRESENT",        # CPBSA
        "",                                                   # CPBSACL → no CLMS dataset
        "EO:CLMS:DAT:CPSCT_10M_YEARLY_2017_PRESENT",        # CPSCT
        "EO:CLMS:DAT:CPSCE_10M_YEARLY_2017_PRESENT",        # CPSCE
        "EO:CLMS:DAT:CPSCD_10M_YEARLY_2017_PRESENT",        # CPSCD
        "",                                                   # CPSCDCL → no CLMS dataset
        "EO:CLMS:DAT:CPFLP_10M_YEARLY_2017_PRESENT",        # CPFLP
        "",                                                   # CPFLPCL → no CLMS dataset
        # CPFLD is commented out in hrl_dict — skip
        # "EO:CLMS:DAT:CPFLD_10M_YEARLY_2017_PRESENT",
        "",                                                   # CPFLDCL → no CLMS dataset
        "EO:CLMS:DAT:CPCSY_10M_YEARLY_2017_PRESENT",        # CPCSY
        "EO:CLMS:DAT:CPCST_10M_YEARLY_2017_PRESENT",        # CPCST
    ],
    # --- Tree cover / Forest parent (EO:EEA:DAT:HRL:TCF) ---
    "EO:EEA:DAT:HRL:TCF" => [
        "EO:CLMS:DAT:TCD_10M_YEARLY_2018_PRESENT",          # TCD
        "",                                                   # TCDCL → no CLMS dataset
        "EO:CLMS:DAT:DLT_10M_YEARLY_2018_PRESENT",          # DLT
        "",                                                   # DLTCL → no CLMS dataset
    ],
    # --- Grassland parent (EO:EEA:DAT:HRL:GRA) ---
    "EO:EEA:DAT:HRL:GRA" => [
        "EO:CLMS:DAT:GRA_10M_YEARLY_2017_PRESENT",          # GRA
        "",                                                   # GRACL → no CLMS dataset
        "EO:CLMS:DAT:HER_10M_YEARLY_2017_PRESENT",          # HER
        "EO:CLMS:DAT:PLOUGH_10M_YEARLY_2017_PRESENT",       # PLOUGH
        "EO:CLMS:DAT:GRAME_10M_YEARLY_2017_PRESENT",        # GRAME
        "",                                                   # GRAMECL → no CLMS dataset
        "EO:CLMS:DAT:GRAMD_10M_YEARLY_2017_PRESENT_4_EVENTS", # GRAMD1-4 (single dataset)
    ],
    # --- Confidence layers (no CLMS dataset — return empty) ---
    "CTYCL" => String[],
    "CPMCHCL" => String[],
    "CPMCECL" => String[],
    "CPMCDCL" => String[],
    "CPBSBCL" => String[],
    "CPBSACL" => String[],
    "CPSCDCL" => String[],
    "CPFLPCL" => String[],
    "CPFLDCL" => String[],
    "TCDCL" => String[],
    "DLTCL" => String[],
    "GRACL" => String[],
    "GRAMECL" => String[],
)

"""
    resolve_hrl_dataset_id(id::String)

Resolve a dataset identifier to the corresponding new CLMS dataset IDs.

Returns a `Vector{String}` of new dataset IDs:
- Old parent dataset IDs (e.g. "EO:EEA:DAT:HRL:CRL") → all associated CLMS datasets
- Old product references (e.g. "CTY", "CPMCH") → single-element vector
- New CLMS dataset IDs → single-element vector (passthrough)
- Confidence layer codes with no CLMS equivalent → empty vector

To get confidence layers, query the corresponding base dataset with the
`datasetIdentifier` query parameter.
"""
function resolve_hrl_dataset_id(id::String)::Vector{String}
    # Direct lookup in hrl_dataset_map (old parent IDs and explicit product codes)
    if haskey(hrl_dataset_map, id)
        return hrl_dataset_map[id]
    end
    # Already a new-style CLMS ID
    startswith(id, "EO:CLMS:DAT:") && return [id]
    # Special case: DLT and TCD use 2018, everything else uses 2017
    if id in ("DLT", "TCD")
        return ["EO:CLMS:DAT:$(id)_10M_YEARLY_2018_PRESENT"]
    end
    # Try to match a key from hrl_dict (new-style product code)
    if haskey(hrl_dict, id)
        return ["EO:CLMS:DAT:$(id)_10M_YEARLY_2017_PRESENT"]
    end
    # Try case-insensitive match on hrl_dict keys
    for k in keys(hrl_dict)
        if lowercase(k) == lowercase(id)
            if k in ("DLT", "TCD")
                return ["EO:CLMS:DAT:$(k)_10M_YEARLY_2018_PRESENT"]
            end
            return ["EO:CLMS:DAT:$(k)_10M_YEARLY_2017_PRESENT"]
        end
    end
    error("Cannot resolve HRL dataset id: '$id'. Must be an old HRL ID (e.g. \"EO:EEA:DAT:HRL:CRL\"), a new CLMS ID (e.g. \"EO:CLMS:DAT:CTY_10M_YEARLY_2017_PRESENT\"), or a product code from HRL.hrl_dict (e.g. \"CTY\", \"CPCSY\", \"CPMCH\").")
end

# Internal helper: download from a vector of already-resolved CLMS dataset IDs
# (skips the resolve step to avoid double-resolution)
function _hrl_download_resolved(
    resolved_ids::Vector{String},
    coord::Tuple,
    yr::Union{String,Int};
    extent=0.05,
    kwargs...,
)::Vector{String}

    hda_client = hrl_hda()
    json = pyimport("json")
    (lon, lat) = coord
    fnames = Vector{String}()

    for ds_id in resolved_ids
        isempty(ds_id) && continue  # skip confidence layers

        query = """
        {
          "dataset_id": "$(ds_id)",
          "bbox": [
            $(lon-extent),
            $(lat-extent),
            $(lon+extent),
            $(lat+extent)
          ],
          "startdate": "$(DateTime(yr)).000Z",
          "enddate": "$(DateTime(yr) + Day(1)).000Z",
          "itemsPerPage": 200,
          "startIndex": 0
        }
        """
        matches = hda_client.search(json.loads(query))
        try
            for idx in pylist(0:(length(matches)-1))
                item = matches.results[idx]
                filename = pyconvert(String, item["id"]) * ".zip"
                zip_path = joinpath(hrl_localpath, filename)
                if !isfile(zip_path)
                    matches[idx].download(download_dir=hrl_localpath)
                end
                push!(fnames, filename)
            end
        catch
            @warn "No matches found for $(ds_id) at $(coord) in year $(yr)"
        end
    end
    return fnames
end


"""
`hrl_hda()`

Sets up client to Harmonised Data Access API on WEkEO
Returns HDA client as Python Object
"""
function hrl_hda()
    # Harmonised Data Access API on WEkEO
    hda = try
        pyimport("hda") #("Client", "Configuration")
    catch
        CondaPkg.add("numpy")
        CondaPkg.add("hda")
        pyimport("hda")
    end
    # Configure user's credentials with a .hdarc file
    hda_client = try
        hda.Client()
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
function hda_datasets(; filteron="EO:CLMS:DAT")
    hda_client = hrl_hda()
    datasets = hda_client.datasets()
    dataset_ids = [item["dataset_id"] for item in datasets]
    # convert to julia
    pyconvert(Vector{String}, dataset_ids) |> filter(x -> occursin(filteron, x))
end

"""
`hda_products(dataset_id)`

Lists all valid productType for `dataset_id` on WEkEO.
Returns a vector of strings.

For new CLMS datasets (EO:CLMS:DAT:*) each dataset has a single product,
so this returns ["default"]. For old HRL datasets the full list is returned.
"""
function hda_products(dataset_id)
    hda_client = hrl_hda()
    metadata = hda_client.metadata(dataset_id)
    props = metadata["properties"]["productType"]
    if !haskey(props, "oneOf")
        # New CLMS datasets have a single product (no oneOf enum)
        return ["default"]
    end
    map(
        x -> x["const"],
        pyconvert(Vector{Dict}, props["oneOf"]),
    )
end

"""
`hrl_info(dataset_id)`

Returns full metadata for a dataset (resolves old HRL IDs to new CLMS IDs first).
When the old ID maps to multiple CLMS datasets, returns metadata for the primary one.
Returns a Dict of the dataset metadata.
"""
function hrl_info(dataset_id)
    resolved = resolve_hrl_dataset_id(dataset_id)
    ds_id = first(resolved)
    hda_client = hrl_hda()
    metadata = hda_client.dataset(ds_id)
    pyconvert(Dict, metadata)
end

# datasetIdentifier (legacy alias)
const hda_info = hrl_info


"""
`hda_years(dataset_id)`

Lists all valid years or epochs for `dataset_id` on WEkEO.
Returns a vector of strings.

For new CLMS datasets (EO:CLMS:DAT:*) the year range is encoded in the
dataset ID itself (e.g. "CTY_10M_YEARLY_2017_PRESENT"), so this returns
["2017_PRESENT"] (extracted from the dataset ID). For old HRL datasets
the full list from the API is returned.
"""
function hda_years(dataset_id)
    hda_client = hrl_hda()
    metadata = hda_client.metadata(dataset_id)
    props = metadata["properties"]["year"]
    if !haskey(props, "oneOf")
        # New CLMS datasets: extract year range from dataset ID
        if occursin("2017_PRESENT", dataset_id)
            return ["2017_PRESENT"]
        elseif occursin("2018_PRESENT", dataset_id)
            return ["2018_PRESENT"]
        elseif occursin("2019_PRESENT", dataset_id)
            return ["2019_PRESENT"]
        elseif occursin("2020_PRESENT", dataset_id)
            return ["2020_PRESENT"]
        end
        return ["PRESENT"]
    end
    map(x -> x["const"], pyconvert(Vector{Dict}, props["oneOf"]))
end


"""
`hrl_download(dataset_id::String, coordinates::Tuple, yr::Union{String,Int}; extent=0.05, product::Union{Nothing,String}=nothing)`
`hrl_download(dataset_id::Vector{String}, coordinates::Tuple, yr::Union{String,Int}; extent=0.05, product::Union{Nothing,String}=nothing)`

Downloads Copernicus Land Monitoring Service High Resolution Layer products from WEkEO using the Harmonised Data Access API.

`dataset_id` can be:
- Old HRL dataset IDs (e.g. "EO:EEA:DAT:HRL:CRL") — resolved to ALL associated CLMS datasets
- New CLMS dataset IDs (e.g. "EO:CLMS:DAT:CTY_10M_YEARLY_2017_PRESENT")
- Product codes from `hrl_dict` (e.g. "CTY", "CPCSY", "CPSCT")

Searches for all tiles near `coordinates` with (lon, lat) +- `extent`, for year `yr`.
Data are downloaded at `hrl_localpath` only if they do not already exist.

When requesting an old parent dataset (CRL, TCF, GRA), all corresponding CLMS
datasets are queried. Confidence layers (CTYCL, TCDCL, etc.) have no dedicated
CLMS dataset — use the base dataset instead.

Returns a Vector of filenames.

See also `hda_datasets`, `hda_products`, `hda_years`
"""
function hrl_download(
    dataset_id::String,
    coord::Tuple,
    yr::Union{String,Int};
    kwargs...,
)
    resolved = resolve_hrl_dataset_id(dataset_id)
    isempty(resolved) && return String[]
    _hrl_download_resolved(resolved, coord, yr; kwargs...)
end

function hrl_download(
    datasets::Vector{String},
    coord::Tuple,
    yr::Union{String,Int};
    kwargs...,
)
    all_resolved = String[]
    for dataset_id in datasets
        append!(all_resolved, resolve_hrl_dataset_id(dataset_id))
    end
    _hrl_download_resolved(all_resolved, coord, yr; kwargs...)
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
    "GRAME" => "Grassland Mowing Events",
    "GRAMECL" => "Grassland Mowing Events Confidence Layer",
    "GRAMD1" => "Grassland Mowing Dates of first mowing event",
    "GRAMD2" => "Grassland Mowing Dates of second mowing event",
    "GRAMD3" => "Grassland Mowing Dates of third mowing event",
    "GRAMD4" => "Grassland Mowing Dates of fourth mowing event",
)

const hrl_legends = Dict(
    "MCTY" => SortedDict(1 => "No cropland", 2 => "Annual crop", 3 => "Permanent crop"),
    "CTY" => SortedDict(
        0 => "No cropland",
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
        65535 => "outside area",
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
        65535 => "Outside area",
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
    "CPFLP" =>
        SortedDict(0 => "No fallow land", 1 => "Fallow land", 65535 => "Outside area"),
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
        65535 => "Outside area",
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
    "GRAME" => SortedDict(
        0 => "no mowing detected",
        1 => "1 mowing event detected",
        2 => "2 mowing events",
        3 => "3 mowing events",
        4 => "4 mowing events",
        253 => "all non-herbaceous areas",
        255 => "outside area",
    ),
    "GRAMECL" => "0-100: Mowing detection confidence. 253: All non-mowing areas. 255: outside area",
    "GRAMD1" => "0: no mowing detected. 1-366: Start (DOY) of each mowing event. 65533: all non-herbaceous areas. 65535: outside area",
    "GRAMD2" => "0: no mowing detected. 1-366: Start (DOY) of each mowing event. 65533: all non-herbaceous areas. 65535: outside area",
    "GRAMD3" => "0: no mowing detected. 1-366: Start (DOY) of each mowing event. 65533: all non-herbaceous areas. 65535: outside area",
    "GRAMD4" => "0: no mowing detected. 1-366: Start (DOY) of each mowing event. 65533: all non-herbaceous areas. 65535: outside area",
)

"""
Extract the year from an HRL downloaded filename.

Handles both old format ("CTY_2023_tile09VFP.zip" → "2023")
and new format ("CLMS_HRLVLCC_CTY_S2024_R10m_...zip" → "2024").
"""
function _hrl_year_from_fname(fname::String)::String
    # Strip .zip if present
    base = endswith(fname, ".zip") ? fname[1:(end-4)] : fname
    parts = split(base, "_")
    # New format: CLMS_HRLVLCC_CTY_S2024_R10m_...
    # Index 4 is "S2024" — strip the S prefix
    if length(parts) >= 4 && occursin(r"^S\d{4}$", parts[4])
        return parts[4][2:end]
    end
    # Old format: CTY_2023_tile09VFP or similar
    for p in parts
        if occursin(r"^\d{4}$", p)
            return p
        end
    end
    return "unknown"
end

"""
`hrl_warp(event::ARCEMEAnalysis.Event; batch="ARCEME-DC-6", dataset_id::Union{String, Vector{String}}=["EO:EEA:DAT:HRL:CRL", "EO:EEA:DAT:HRL:TCF", "EO:EEA:DAT:HRL:GRA"])`
`hrl_warp(cubename, fnames::Vector{String}; batch="ARCEME-DC-6")`
`hrl_warp(cubename::String; batch="ARCEME-DC-6", dataset_id::Union{String, Vector{String}}=["EO:EEA:DAT:HRL:CRL", "EO:EEA:DAT:HRL:TCF", "EO:EEA:DAT:HRL:GRA"])`

Reprojects and clips files in `fnames` to match `event`'s coordinate reference system and extent,
and saves them as a Zarr group in `batch-HRL` under `ARCEMEAnalysis.local_cubepath`.

Old dataset IDs (e.g. "EO:EEA:DAT:HRL:CRL") are resolved automatically for backward compatibility.
"""
function hrl_warp(cubename, fnames::Vector; batch="ARCEME-DC-6", suffix="", zip=true)
    ds = arceme_open(cubename; batch)
    event_date = arceme_eventdate(ds)
    res = ds.properties["resolution"]
    tmp = tempname(; suffix=".tif")

    outputzarr = joinpath(ARCEMEAnalysis.local_cubepath, "$(batch)-HRL", cubename)

    for hrl in keys(hrl_dict)
        hrl_fnames = fnames[map(contains("$(hrl)_"), fnames)]
        if isempty(hrl_fnames)
            continue
        end
        # use gdal distribution shipped with (Arch)GDAL.jl rather than system One
        try
            run(
                `$(AG.GDAL.gdalwarp_exe()) 
                -t_srs $(AG.toPROJ4(AG.importEPSG(ds.properties["epsg"]))) 
                -tr $res $res
                -te $(minimum(ds.x)) $(minimum(ds.y)) $(maximum(ds.x)+res) $(maximum(ds.y)+res) 
                -r near
                $(["/vsizip/$(joinpath(hrl_localpath, hrl_fnames[i], "$(split(hrl_fnames[i],".")[1]).tif") )" for i in eachindex(hrl_fnames)]) 
                $tmp 
                -overwrite 
                `
            )
        catch
            run(
                `$(AG.GDAL.gdalwarp_exe()) 
                -t_srs $(AG.toPROJ4(AG.importEPSG(ds.properties["epsg"]))) 
                -tr $res $res
                -te $(minimum(ds.x)) $(minimum(ds.y)) $(maximum(ds.x)+res) $(maximum(ds.y)+res) 
                -r near
                $(["/vsizip/$(joinpath(hrl_localpath, hrl_fnames[i], "$(split(hrl_fnames[i],".")[1])","$(split(hrl_fnames[i],".")[1]).tif") )" for i in eachindex(hrl_fnames)]) 
                $tmp 
                -overwrite 
                `
            )
        end
        c = Cube(tmp)
        # write to zarr
        props = Dict(
            "long_name" => hrl_dict["$(split(hrl_fnames[1],"_")[3])"],
            "year_coverage" => _hrl_year_from_fname(hrl_fnames[1]),
            "legend" => hrl_legends["$(split(hrl_fnames[1],"_")[3])"],
            "epsg" => ds.properties["epsg"],
            "resolution" => 10,
            "unit" => "meter",
            "source" => "$(join(hrl_fnames, ", ")) © European Union, Copernicus Land Monitoring Service 2025, European Environment Agency (EEA)",
            "processing" => "Downloaded from wEkEO. Reprojected to EPSG $(ds.properties["epsg"]) using nearest neighour, mosaicked and clipped to AOI with gdalwarp",
            "processed_by" => "M. Weynants, MPI-BGC",
            "processed_date" => string(Date(now())),
        )
        c2 = YAXArray((ds.x, ds.y), c.data, props)
        nt = (Symbol("$(split(hrl_fnames[1],"_")[3])$suffix") => c2,)
        tds = Dataset(; nt...)
        tds = setchunks(tds, (500, 500))
        savedataset(
            tds,
            path=outputzarr,
            append=true,
        )
    end
    # zip ?
    !zip && return outputzarr
    outputzip = "$outputzarr.zip"
    isfile(outputzip) && rm(outputzip)
    ARCEMEAnalysis.zip_dir(outputzip, outputzarr)
    rm(outputzarr, recursive=true)
    return outputzip
end

hrl_warp(event::ARCEMEAnalysis.Event; kwargs...) =
    hrl_warp(arceme_cubename(event); kwargs...)

function hrl_warp(
    cubename::String;
    batch="ARCEME-DC-6",
    dataset_id::Union{String,Vector{String}}=[
        "EO:EEA:DAT:HRL:CRL",
        "EO:EEA:DAT:HRL:TCF",
        "EO:EEA:DAT:HRL:GRA",
    ],
    yr=nothing,
    kwargs...,
)
    ds = arceme_open(cubename; batch)
    coord = arceme_coordinates(ds)
    if isnothing(yr)
        yr = year(arceme_eventdate(ds))
    end
    fnames = hrl_download(dataset_id, coord, yr)
    hrl_warp(cubename, fnames; batch, kwargs...)
end

"""
     arceme_stats(ds::YAXArrays.Dataset, name; classes=nothing)
     arceme_stats(ev::Event, name; batch="ARCEME-DC-8", classes=nothing)

Get the  statistics of layer `name` for the ARCEME data cube `ds`.
Default classes are extracted from `hrl_legends[name]`. Alternative classes can be provided as a dictionary.
"""
function arceme_stats(ds, name; classes=nothing)
    if isnothing(classes)
        classes = hrl_legends[name]
    end
    cdr = counter(ds[name])
    map(collect(classes)) do (k, v)
        count = get(cdr, k, 0)
        (
            key=k,
            class=v,
            count=count,
            fraction=count / (length(ds.x) * length(ds.y)),
        )
    end
end
arceme_stats(ev::ARCEMEAnalysis.Event, name; batch="ARCEME-DC-8", classes=nothing) =
    arceme_stats(arceme_open(ev, batch=batch), name; classes=classes)

function ctykeymap(k)
    k == 0 && return 1
    k > 3000 && return 1
    ctL1 = k ÷ 1000
    ctL2 = (k - 1000 * ctL1) ÷ 100
    ctL3 = (k - 1000 * ctL1 - 100 * ctL2) ÷ 10
    ctL1 == 1 && ctL2 == 1 && return ctL3 + 1 # cereals
    ctL1 == 1 && ctL2 == 2 && return ctL3 + 6 # vegetables and pulses
    ctL1 == 1 && ctL2 == 3 && return ctL3 + 8 # roots (potato, sugar beet)
    ctL1 == 1 && ctL2 == 4 && return ctL3 + 10 # oil (sunflower, soy, rapeseed)
    ctL1 == 2 && ctL2 == 1 && return ctL3 + 15 # fiber (flax, cotton, hemp)
    ctL1 == 2 && ctL2 == 2 && return ctL2 + 14 # grape, olive
    ctL1 == 2 && ctL2 == 3 && return ctL3 + 16 # fruit, nuts
end

function mctykeymap(k)
    k == 0 && return 1 # non crop
    k == 65535 && return 1
    ctL1 = k ÷ 1000
    ctL2 = (k - 1000 * ctL1) ÷ 100
    ctL3 = (k - 1000 * ctL1 - 100 * ctL2) ÷ 10
    ctL1 == 1 && return 2 # annual crops
    # ctL1 == 2 && ctL2 == 1 && return 2 # annual fiber
    ctL1 == 2 && return 3 # permanent crop
    k == 3100 && return 2 # Unclassified annual crop
    k == 3200 && return 3 # Unclassified permanent crop
end


function get_soilcover(event_date, ds)
    amce = ARCEMEAnalysis.CubePlots.get_mce(ds)
    amch = ARCEMEAnalysis.CubePlots.get_mch(ds)
    asce = ARCEMEAnalysis.CubePlots.get_sce(ds)
    ascd = ARCEMEAnalysis.CubePlots.get_val(ds, layer="CPSCD")
    absb = ARCEMEAnalysis.CubePlots.get_val(ds, layer="CPBSB")
    absa = ARCEMEAnalysis.CubePlots.get_val(ds, layer="CPBSA")
    soilcover = broadcast(amce, amch, asce, ds.CPSCT, ascd, absb, absa) do mce, mch, sce, sct, scd, bsb, bsa
        (ismissing(mce) || ismissing(mch)) && return UInt8(255)
        if event_date >= mce
            if event_date <= mch
                return UInt(1)
            else
                ismissing(bsa) && return missing
                if event_date < mch + Day(bsa)
                    return UInt8(0)
                else
                    ismissing(scd) && return missing
                    if sct in (1, 2) && event_date > (sce + Day(scd)) # summer secondary crop
                        return UInt8(0)
                    else
                        return UInt8(1)
                    end
                end
            end
        else # event_date < mce
            ismissing(bsb) && return missing
            if event_date >= (mce-Day(bsb))
                return UInt8(0)
            else
                ismissing(scd) && return missing
                if sct in (3, 4) && event_date < sce  # winter secondary crop
                    return UInt8(0)
                else
                    return UInt8(1)
                end
            end
        end
    end
    return soilcover
end
export hrl_download,
    hrl_set_localpath, hda_datasets, hda_products, hda_years, hrl_warp, arceme_stats, hrl_info, resolve_hrl_dataset_id
end #module