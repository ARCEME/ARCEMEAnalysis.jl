module FingerprintPlots
using Colors: RGBA
using Makie
using Dierckx
using YAXArrays
using Dates
using ..ARCEMEAnalysis
using ..ARCEMEAnalysis: local_cubepath, arceme_legends, Event
using DataStructures: counter, OrderedDict
using Statistics: mean
using GeoMakie
import DimensionalData as DD

function filter_spline(xout, samplets, x, nout)
    missmask = findall(x -> !ismissing(x) && isfinite(x), samplets)
    isempty(missmask) && return xout .= missing
    xgood = x[missmask]
    ygood = samplets[missmask]

    spl = nothing
    nknots = 19

    kn = compute_knots(xgood, nknots)

    spl = try
        Spline1D(xgood, ygood, kn)
    catch
        return xout .= NaN
    end

    xout .= evaluate(spl, range(0, 2, length=nout))
end




pairnow = arceme_validpairs()[4]



function most_common_class(pairnow)
    ds1, ds2 = arceme_open.(pairnow)
    #Find maximum landcover in both cubes
    classfrac1, classfrac2 = (ds1.class_fraction[:], ds2.class_fraction[:])
    (_, i1), (_, i2) = findmax.((classfrac1, classfrac2))
    if classfrac1[i1] * classfrac2[i1] > classfrac1[i2] * classfrac2[i2]
        ds1.class.val[i1]
    else
        ds1.class.val[i2]
    end
end

function pairwise_without_missings(x, y)
    out = Point2f[]
    for i in eachindex(x, y)
        if !ismissing(y[i]) && isfinite(y[i])
            push!(out, Point2f(x[i], y[i]))
        end
    end
    out
end

function open_plot_data(ev, class; strata="ESA_LC")
    @info "Opening $ev"
    ds = arceme_open(ev)
    firstday = arceme_eventdate(ev) - Year(1)
    lastday = arceme_eventdate(ev) + Year(1)
    nms = Millisecond(lastday - firstday).value
    linscal = d -> Millisecond(d - firstday).value * 2 / nms
    x_s2 = map(linscal, ds.time_sentinel_2_l2a.val)
    x_s1 = map(linscal, ds.time_sentinel_1_rtc.val)
    outax = YAXArrays.time(range(firstday, lastday, step=Day(1)))

    s2 = ds.s2_indices[class=DD.At(class)]
    s1 = ds.s1_indices[class=DD.At(class)]

    r_s2 = xmap(filter_spline, s2 ⊘ :time_sentinel_2_l2a, output=XOutput(outax), function_args=(x_s2, length(outax),))
    r_s1 = xmap(filter_spline, s1 ⊘ :time_sentinel_1_rtc, output=XOutput(outax), function_args=(x_s1, length(outax),))

    dataplots = OrderedDict()
    interplots = OrderedDict()

    for band in ds.band_s2

        interplots[band] = pairwise_without_missings(range(0, 2, length=length(outax)), r_s2[band_s2=DD.At(band)].data[:])
        dataplots[band] = pairwise_without_missings(x_s2, s2[band_s2=DD.At(band)].data[:])
    end

    for band in ds.band_s1
        interplots[band] = pairwise_without_missings(range(0, 2, length=length(outax)), r_s1[band_s1=DD.At(band)].data[:])
        dataplots[band] = pairwise_without_missings(x_s1, s1[band_s1=DD.At(band)].data[:])
    end

    classfracs = ds.class_fraction.data[:]
    sfracs = sortperm(classfracs, rev=true)
    strs = map(sfracs[1:3]) do ifrac
        string(collect(values(arceme_legends[strata]))[ifrac], " => ", round(lcfracs[ifrac] * 100), "%")
    end
    classdesc = join(strs, ", ")


    (dataplots, interplots, classdesc)
end



function fingerprint_plot(pairid=1; class=nothing, plotdict=IdDict(), interactive=true, strata="ESA_LC")

    ipair = Observable{Int}(pairid)
    alleventpairs = arceme_validpairs()

    if class === nothing
        class = most_common_class(alleventpairs[pairid])
    end

    class_choos = Observable{Any}(class)

    d = lift(ipair, class_choos) do i, class
        ev = alleventpairs[i][1]
        open_plot_data(ev, class)
    end
    dhp = lift(ipair, class_choos) do i, class
        ev = alleventpairs[i][2]
        open_plot_data(ev, class)
    end
    fig = Figure(size=(1600, 2600))

    allbands = collect(keys(d[][1]))

    allax = Dict()

    for (i_d, o) in enumerate((d, dhp))
        colcur = (:red, :blue)[i_d]
        dataplots = lift(i -> i[1], o)
        interplots = lift(i -> i[2], o)
        i_pl = 1
        for index_choos in allbands
            ax = Axis(fig[i_pl, i_d], title=index_choos)
            allax[(i_pl, i_d)] = ax
            data_interp = lift(i -> i[index_choos], interplots)
            data_orig = lift(i -> i[index_choos], dataplots)
            lines!(ax, data_interp, color=colcur)
            p = plot!(ax, data_orig, color=colcur)
            plotdict[p] = (index_choos, i_d)
            i_pl += 1
        end
    end

    npl = maximum(first, keys(allax))
    for ipl in 1:npl
        linkyaxes!(allax[(ipl, 1)], allax[(ipl, 2)])
    end
    ov1 = Axis(fig[5:6, 3], aspect=DataAspect(), title="Drought")
    img1 = lift((y, i) -> arceme_representative_image(alleventpairs[i][1], class=y), class_choos, ipair)
    heatmap!(ov1, img1)
    ov2 = Axis(fig[7:8, 3], aspect=DataAspect(), title="Drought + HP")
    img2 = lift((y, i) -> arceme_representative_image(alleventpairs[i][2], class=y), class_choos, ipair)
    heatmap!(ov2, img2)


    xpos = lift(i -> [alleventpairs[i][1].longitude, alleventpairs[i][2].longitude], ipair)
    ypos = lift(i -> [alleventpairs[i][1].latitude, alleventpairs[i][2].latitude], ipair)
    img = YAXArray((lon(-179.75:0.5:179.75), lat(-89.75:0.5:89.75)), rotr90(GeoMakie.earth()))
    gax = GeoAxis(fig[3:4, 3], source="+proj=longlat", dest="+proj=longlat")
    heatmap!(gax, img; interpolate=false)
    scatter!(gax, xpos, ypos, color=[:red, :blue])

    # xlims!(gax, lift(xpos -> minimum(xpos) - 20, xpos), lift(xpos -> maximum(xpos) + 20, xpos))
    # ylims!(gax, lift(ypos -> minimum(ypos) - 20, ypos), lift(ypos -> maximum(ypos) + 20, ypos))
    if interactive
        menupair = Textbox(fig, placeholder=string(ipair[]), validator=Int, tellwidth=false)
        #menupair = Slider(fig, range=1:length(alleventpairs), startvalue=pairid, update_while_dragging=false)

        menuclass = Menu(fig, options=collect(values(arceme_legends[strata])), default=class)
        fig[1:2, 3] = vgrid!(
            menupair,
            Label(fig, lift(i -> i[3], d), color=:red),
            Label(fig, lift(i -> i[3], dhp), color=:blue),
            menuclass,)
        on(menuclass.selection) do s
            class_choos[] = s
        end
        on(menupair.stored_string) do v
            ipair[] = parse(Int, v)
            for ipl in 1:npl
                autolimits!(allax[(ipl, 1)])
            end
            autolimits!(ov1)
            autolimits!(ov2)
        end
    end
    fig
end


arceme_representative_image(ev::Event; index_use="NDVI", brighten_factor=2, class=nothing) =
    arceme_representative_image(arceme_open(ev); index_use, brighten_factor, class)

"""
    arceme_representative_image(ev)

Creates an RGB image highlighting on a cloud-free time step with high NDVI for a certain land cover type 
"""

function arceme_representative_image(ds; index_use="NDVI", brighten_factor=2, class=nothing, strata="ESA_LC")
    if class === nothing
        classfrac = ds.class_fraction.data[:]
        _, iclass = findmax(classfrac)
        class = ds.class.val[iclass]
    end
    ndvi = ds.s2_indices[band=DD.At(index_use), class=DD.At(class)]
    allndviranks = sortperm(ndvi.data[:], rev=true)
    repr_ind = 0
    cloudfrac = ds.cloud_fraction.data[:]
    for i in 1:length(allndviranks)
        icandidate = findfirst(==(i), allndviranks)
        if !ismissing(ndvi[icandidate]) && cloudfrac[icandidate] < 0.05
            repr_ind = icandidate
            break
        end
    end
    if repr_ind == 0
        return YAXArray((ds.x, ds.y), fill(RGBA(0.0, 0.0, 0.0, 1.0), 1000, 1000), Dict{Any,Any}())
    end
    brfilt(a, factor, alpha) = RGBA(clamp(a.r * factor, 0, 1), clamp(a.g * factor, 0, 1), clamp(a.b * factor, 0, 1), alpha)
    iclass = findfirst(==(class), collect(values(arceme_legends[strata])))
    r = broadcast(arceme_rgb(ds)[time_sentinel_2_l2a=repr_ind], ds.ESA_LC[:, :, 1], brighten_factor) do col, class, f
        classv = ARCEMEAnalysis.lckeymap(class)
        alpha = classv == iclass ? 1.0 : 0.5
        brfilt(col, f, alpha)
    end
    r[:, :]
end

function compute_knots(xgood, ntarget=19)
    xtarget = range(0.0, 2.0, length=ntarget)
    kn = filter(x -> first(xgood) < x < last(xgood), xtarget)
    knotstoreplace = []
    for i in 1:(length(kn)-1)
        ix = searchsortedfirst(xgood, kn[i])
        ix2 = searchsortedfirst(xgood, kn[i+1])
        if ix2 - ix <= 1
            if !isempty(knotstoreplace) && last(last(knotstoreplace)) == i
                push!(last(knotstoreplace), i + 1)
            else
                push!(knotstoreplace, [i, i + 1])
            end
        end
    end
    if isempty(knotstoreplace)
        return kn
    else
        newknots = map(knotstoreplace) do group
            mean(view(kn, group))
        end
        return sort([kn[setdiff(1:length(kn), reduce(vcat, knotstoreplace))]; newknots])
    end
end

export fingerprint_plot
end