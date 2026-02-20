import Dates: Millisecond
import YAXArrays: YAXArray
"""
    arceme_s1_position(ds)

Analyses the sentinel 1 time stamps of the datasets and groups them into series of time stamps whose difference
is always a multiple of 6 days, so we can assume that all images with the same group tag are retrieved from the 
same position. Store the groups in the dataset with the name `s1_postion`.
"""
function arceme_s1_position(ds)
    ts = ds.time_sentinel_1_rtc.val
    sixdays = 6*24*60*60*1000
    group_offsets = [0]
    groups = [[1]]
    for its in 2:length(ts)
        groupfound = false
        for igroup in eachindex(groups)
            fac = (Millisecond(ts[its]-ts[1]).value-group_offsets[igroup])/sixdays
            if abs(round(fac)-fac) < 1e-6
                push!(groups[igroup],its)
                groupfound=true
                break
            end
        end
        if !groupfound
            newoffset = Dates.Millisecond(ts[its]-ts[1]).value-floor(Int,fac)*sixdays
            push!(group_offsets,newoffset)
            push!(groups,[its])
        end
    end
    positions = zeros(Int,length(ts))
    for i in 1:length(groups)
        positions[groups[i]] .= i
    end
    positions
    ds.cubes[:s1_position] = YAXArray((ds.time_sentinel_1_rtc,),positions)
    nothing
end

_to_db(x) = 10.0 * log10(x)

function arceme_radar_db(ds)
    ds.cubes[:vv_db] = _to_db.(ds.vv)
    ds.cubes[:vh_db] = _to_db.(ds.vh)
end