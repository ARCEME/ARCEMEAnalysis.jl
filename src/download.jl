using StringViews
using ProgressMeter
using HTTP, JSON
using ZipArchives
using Preferences: @set_preferences!, @load_preference
export arceme_set_localpath

"""
    arceme_set_localpath(path)

Sets the path to a local directory where downloaded copies of the arceme cubes are stored. 
"""
function arceme_set_localpath(path)
    @set_preferences!("arceme_localpath" => path)
end

const local_cubepath = @load_preference("arceme_localpath")

function arceme_download_batch(batch="ARCEME-DC-6")
    if local_cubepath === nothing
        error("You need to set the local cube path first. Please run `arceme_localpath(path)` first.")
    end
    aresp = HTTP.get("https://s3.waw3-2.cloudferro.com/swift/v1/$batch/",query=Dict("format"=>"json","delimiter"=>"/"))
    allarrays = map(i->strip(i["subdir"],'/'),JSON.parse(aresp.body))

    @showprogress for current_cube in allarrays
        lresp = HTTP.get("https://s3.waw3-2.cloudferro.com/swift/v1/$batch/",query=Dict("prefix"=>current_cube))
        files_in_cube = split(StringView(lresp.body),"\n")

        outfilename = joinpath(local_cubepath,batch,string(current_cube,".zip"))
        ZipWriter(outfilename) do w
            for f in files_in_cube
                resp = HTTP.get("https://s3.waw3-2.cloudferro.com/swift/v1/$batch/$f");
                f2 = joinpath(splitpath(f)[2:end]...)
                zip_newfile(w, f2)
                write(w,resp.body)
            end
        end
    end
end