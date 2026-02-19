module HRL
import CondaPkg
using  PythonCall: pyimport, pyconvert, pylist
import GeoJSON
using Preferences: @set_preferences!, @load_preference

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


dataset_id = "EO:EEA:DAT:HRL:CRL";yr="2023";coord=(-1.0, 38.6);product="Crop Types";extent=0.05
"""
`hrl_download(dataset_id::String, coordinates::Tuple, yr::Union{String,Int}; extent=0.05, hrl_localpath=hrl_localpath, product::{Union{Nothing,String}}=nothing)`

Downloads Copernicus Land Monitoring Service High Resolution Layer products from WEkEO using the Harmonised Data Access API. 
  `dataset_id` must be a valid dataset (e.g., Cropland: "EO:EEA:DAT:HRL:CRL"; Tree cover and Forest: "EO:EEA:DAT:HRL:TCF"; 
  Grassland: "EO:EEA:DAT:HRL:GRA").
  Search for all tiles near `coordinates` with (lon, lat) +- `extent`, for year `yr`. Some products are only available for three-years epochs.
  Data are downloaded at `hrl_localpath` only if they do not already exist.
  By default, all products from `dataset_id` are returned. To search for a specific product, set `product` to a valid productType (e.g. "Crop Types")
  
Returns a Vector of filenames.

See also `hda_datasets`, `hda_products`, `hda_years`
"""
function hrl_download(dataset_id::String, coord::Tuple, yr::Union{String,Int}; extent=0.05, hrl_localpath=hrl_localpath, product::Union{Nothing,String}=nothing)
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
    for i in pylist(findall(todownload).-1)
      println("Downloading $(matches[i]["id"])")
      matches[i].download(download_dir=hrl_localpath) 
    end
    return fnames
end

export hrl_download, hrl_set_localpath, hda_datasets, hda_products, hda_years
end #module