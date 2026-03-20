module FingerprintPlots
using Colors:RGBA
using Makie
using Dierckx
using YAXArrays
using Dates
using ..ARCEMEAnalysis
using ..ARCEMEAnalysis: local_cubepath
using DataStructures: counter
using Statistics: mean
using GeoMakie
import DimensionalData as DD

function filter_spline(xout, samplets,x,nout)
    missmask = findall(x->!ismissing(x) && isfinite(x),samplets)
    isempty(missmask) && return xout .= missing
    xgood = x[missmask]
    ygood = samplets[missmask]
    
    spl = nothing
    nknots = 19

    kn = compute_knots(xgood,nknots)

    spl = try
        Spline1D(xgood,ygood,kn)
    catch
        return xout .= NaN
    end

    xout .= evaluate(spl,range(0,2,length=nout))
end


function arceme_representative_image(ev; index_use="kNDVI", brighten_factor=2, lc = nothing)
    ds = arceme_open(ev)
    fp = arceme_open_fingerprint(ev)
    sumcl = xmap(ds.cloud_mask ⊘ (:x,:y),ds.SCL ⊘ (:x,:y),inplace=false,output=XOutput(outtype=Float64)) do cl,scl
        sum(i->ARCEMEAnalysis._is_cloud(i...),zip(cl,scl))/length(cl)
    end
    sumcl = sumcl[1,1,:]
    alllc = if lc === nothing 
        xmap(fp.s2_indices[band=DD.At(index_use)] ⊘ :lc, inplace=false) do vals
            all(ismissing,vals) && return NaN
            mean(skipmissing(vals))
        end
    else
        fp.s2_indices[band=DD.At(index_use),lc=DD.At(lc)]
    end
    alllcranks = sortperm(alllc.data[:],rev=true)
    repr_ind = 0
    for i in 1:length(alllcranks)
        icandidate = findfirst(==(i),alllcranks)
        if sumcl.data[icandidate]< 0.05
            repr_ind = icandidate
            break
        end
    end
    brfilt(a,factor,alpha) = RGBA(clamp(a.r*factor,0,1),clamp(a.g*factor,0,1),clamp(a.b*factor,0,1),alpha)
    ilc = findfirst(==(lc),collect(values(arceme_classes)))
    r = broadcast(arceme_rgb(ds)[time_sentinel_2_l2a=repr_ind],ds.ESA_LC[:,:,1],brighten_factor) do col,lc,f
        lcv = ARCEMEAnalysis.lckeymap(lc)
        alpha = lcv == ilc ? 1.0 : 0.5
        brfilt(col,f,alpha)
    end
    r[:,:]
end

pairnow = arceme_validpairs()[4]



function most_common_lc(pairnow)
    c1,c2 = map(pairnow) do ev
        ds = arceme_open(ev)
        counter(ds.ESA_LC.data[:,:,1])
    end
    #Find maximum landcover in both cubes
    (_,i1),(_,i2) = findmax.((c1,c2))
    if c1[i1] + c2[i1] > c1[i2] + c2[i2]
        arceme_classes[i1]
    else
        arceme_classes[i2]
    end
end



function fingerprint_plot(pairnow; lc=most_common_lc(pairnow))
    lc_choos = lc
    d,dhp = map(pairnow) do ev
        ds = arceme_open_fingerprint(ev)
        firstday = arceme_eventdate(ev) - Year(1)
        lastday = arceme_eventdate(ev) + Year(1)
        nms = Millisecond(lastday-firstday).value
        linscal = d->Millisecond(d-firstday).value*2/nms
        x_s2 = map(linscal,ds.time_sentinel_2_l2a.val)
        x_s1 = map(linscal,ds.time_sentinel_1_rtc.val)
        outax = YAXArrays.time(range(firstday,lastday,step=Day(1)))

        r_s2 = xmap(filter_spline,ds.s2_indices ⊘ :time_sentinel_2_l2a, output=XOutput(outax), function_args=(x_s2, length(outax),))
        r_s1 = xmap(filter_spline,ds.s1_indices ⊘ :time_sentinel_1_rtc, output=XOutput(outax), function_args=(x_s1, length(outax),))

        ds, r_s1, r_s2, x_s1, x_s2
    end

    fig = Figure(size=(1600,2600),title=lc_choos);

    allax = Dict()

    for (i_d,(ds,r_s1,r_s2,x_s1,x_s2)) in enumerate((d,dhp))
        colcur = (:red,:blue)[i_d]
        x,r = x_s2,r_s2
        for (i,index_choos) in enumerate(ds.band_s2.val)
            orig = ds.s2_indices[lc=DD.At(lc_choos),band_s2=DD.At(index_choos)].data[:]
            x_orig = x[(!ismissing).(orig)]
            ax = Axis(fig[i,i_d],title=index_choos)
            allax[(i,i_d)] = ax 
            lines!(ax,range(0,2,length(r.time)),r[lc=DD.At(lc_choos),band_s2=DD.At(index_choos)].data[:],color=colcur);
            plot!(ax,x_orig,orig[(!ismissing).(orig)],color=colcur)
        end
        x,r = x_s1,r_s1

        for (i,index_choos) in enumerate(ds.band_s1.val)
            orig = ds.s1_indices[lc=DD.At(lc_choos),band_s1=DD.At(index_choos)].data[:]
            x_orig = x[(!ismissing).(orig)]
            ax = Axis(fig[i+length(ds.band_s2),i_d],title=index_choos) 
            allax[(i+length(ds.band_s2),i_d)] = ax
            lines!(ax,range(0,2,length(r.time)),r[lc=DD.At(lc_choos),band_s1=DD.At(index_choos)].data[:],color=colcur);
            plot!(ax,x_orig,orig[(!ismissing).(orig)],color=colcur)
        end
    end

    npl = maximum(first,keys(allax))
    for ipl in 1:npl
        linkyaxes!(allax[(ipl,1)],allax[(ipl,2)])
    end
    ax = Axis(fig[4:7,3],aspect=DataAspect(),title="Drought")
    heatmap!(ax,arceme_representative_image(pairnow[1],lc=lc_choos))
    ax = Axis(fig[8:11,3],aspect=DataAspect(),title="Drought + HP")
    heatmap!(ax,arceme_representative_image(pairnow[2],lc=lc_choos))

    x = [pairnow[1].longitude, pairnow[2].longitude]
    y = [pairnow[1].latitude, pairnow[2].latitude]
    img = YAXArray((lon(-179.75:0.5:179.75),lat(-89.75:0.5:89.75)),rotr90(GeoMakie.earth()));
    gax = GeoAxis(fig[1:3,3],source="+proj=longlat",dest = "+proj=longlat")
    heatmap!(gax, img; interpolate = false)
    scatter!(gax,x,y,color=[:red,:blue])
    xlims!(gax,minimum(x)-20,maximum(x)+20)
    ylims!(gax,minimum(y)-20,maximum(y)+20)
    fig
end

function compute_knots(xgood,ntarget=19)
    xtarget = range(0.0,2.0,length=ntarget)
    kn = filter(x->first(xgood) < x < last(xgood), xtarget)
    knotstoreplace = []
    for i in 1:(length(kn)-1)
        ix = searchsortedfirst(xgood,kn[i])
        ix2 = searchsortedfirst(xgood,kn[i+1])
        if ix2-ix <= 1
            if !isempty(knotstoreplace) && last(last(knotstoreplace)) == i
                push!(last(knotstoreplace), i+1)
            else
                push!(knotstoreplace,[i,i+1])
            end
        end
    end
    if isempty(knotstoreplace)
        return kn
    else
        newknots = map(knotstoreplace) do group
            mean(view(kn,group))
        end
        return sort([kn[setdiff(1:length(kn),reduce(vcat,knotstoreplace))];newknots])
    end
end

export fingerprint_plot
end