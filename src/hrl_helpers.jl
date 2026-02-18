module HRL
import CondaPkg
using  PythonCall: pyimport, pyconvert, pylist
import GeoJSON
using Preferences: @set_preferences!, @load_preference
export hrl_download, hrl_set_localpath

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


hrl = [
  (dataset_id = "EO:EEA:DAT:HRL:CRL",),
  (dataset_id = "EO:EEA:DAT:HRL:TCF",),
  (dataset_id = "EO:EEA:DAT:HRL:GRA",),
]



function hrl_download(dataset_id::String, coordinates::Tuple, yr::Union{String,Int}; extent=0.05, hrl_localpath=hrl_localpath)
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
    json = pyimport("json")

    (lon,lat) = coordinates
    fnames = []

    # Build query as a JSON object
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
      "startIndex": 0
    }
    """
    matches = hda_client.search(json.loads(query))
    try
    for item in matches.results
      println(item["id"])
      filename = pyconvert(String, item["id"])*".zip"
      # store info:
      push!(fnames, filename)
    end
    catch
      @warn "No matches found for $(dataset_id) at $(coordinates) in year $(yr)"
      return fnames
    end

    # Download
    @show todownload = [!isfile(joinpath(hrl_localpath, pyconvert(String, item["id"])*".zip")) for  item in matches.results]
    for i in pylist(findall(todownload).-1)
      println("Downloading $(matches[i]["id"])")
      matches[i].download(download_dir=hrl_localpath) 
    end
    return fnames
end

end #module