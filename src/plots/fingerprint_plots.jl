module FingerprintPlots
using Colors: RGBA
using Makie
using Dierckx
using YAXArrays
using Dates
using ..ARCEMEAnalysis
using ..ARCEMEAnalysis: local_cubepath, arceme_classes
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



function most_common_lc(pairnow)
    c1, c2 = map(pairnow) do ev
        ds = arceme_open(ev)
        counter(ds.ESA_LC.data[:, :, 1])
    end
    #Find maximum landcover in both cubes
    (_, i1), (_, i2) = findmax.((c1, c2))
    if c1[i1] + c2[i1] > c1[i2] + c2[i2]
        arceme_classes[i1]
    else
        arceme_classes[i2]
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

function open_plot_data(ev, lc)
    @info "Opening $ev"
    ds = arceme_open_fingerprint(ev)
    firstday = arceme_eventdate(ev) - Year(1)
    lastday = arceme_eventdate(ev) + Year(1)
    nms = Millisecond(lastday - firstday).value
    linscal = d -> Millisecond(d - firstday).value * 2 / nms
    x_s2 = map(linscal, ds.time_sentinel_2_l2a.val)
    x_s1 = map(linscal, ds.time_sentinel_1_rtc.val)
    outax = YAXArrays.time(range(firstday, lastday, step=Day(1)))

    s2 = ds.s2_indices[lc=DD.At(lc)]
    s1 = ds.s1_indices[lc=DD.At(lc)]

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

    (dataplots, interplots)
end


function fingerprint_plot(pairid=1; lc=nothing, plotdict=IdDict())

    ipair = Observable{Int}(pairid)
    alleventpairs = arceme_validpairs()

    if lc === nothing
        lc = most_common_lc(alleventpairs[pairid])
    end

    lc_choos = Observable{Any}(lc)

    d = lift(ipair, lc_choos) do i, lc
        ev = alleventpairs[i][1]
        open_plot_data(ev, lc)
    end
    dhp = lift(ipair, lc_choos) do i, lc
        ev = alleventpairs[i][2]
        open_plot_data(ev, lc)
    end
    fig = Figure(size=(1600, 2600))

    allbands = collect(keys(d[][1]))

    allax = Dict()

    for (i_d, o) in enumerate((d, dhp))
        colcur = (:red, :blue)[i_d]
        dataplots = lift(first, o)
        interplots = lift(last, o)
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
    ax = Axis(fig[5:6, 3], aspect=DataAspect(), title="Drought")
    img1 = lift((y, i) -> arceme_representative_image(alleventpairs[i][1], lc=y), lc_choos, ipair)
    heatmap!(ax, img1)
    ax = Axis(fig[7:8, 3], aspect=DataAspect(), title="Drought + HP")
    img2 = lift((y, i) -> arceme_representative_image(alleventpairs[i][2], lc=y), lc_choos, ipair)
    heatmap!(ax, img2)


    xpos = lift(i -> [alleventpairs[i][1].longitude, alleventpairs[i][2].longitude], ipair)
    ypos = lift(i -> [alleventpairs[i][1].latitude, alleventpairs[i][2].latitude], ipair)
    img = YAXArray((lon(-179.75:0.5:179.75), lat(-89.75:0.5:89.75)), rotr90(GeoMakie.earth()))
    gax = GeoAxis(fig[3:4, 3], source="+proj=longlat", dest="+proj=longlat")
    heatmap!(gax, img; interpolate=false)
    scatter!(gax, xpos, ypos, color=[:red, :blue])
    # xlims!(gax, lift(xpos -> minimum(xpos) - 20, xpos), lift(xpos -> maximum(xpos) + 20, xpos))
    # ylims!(gax, lift(ypos -> minimum(ypos) - 20, ypos), lift(ypos -> maximum(ypos) + 20, ypos))
    menupair = Slider(fig, range=1:length(alleventpairs), startvalue=40, update_while_dragging=false)
    menulc = Menu(fig, options=collect(values(arceme_classes)), default=lc)
    fig[1:2, 3] = vgrid!(
        Label(fig, "Event Pair", width=nothing),
        menupair,
        Label(fig, "Land Cover", width=nothing),
        menulc,
    )
    on(menulc.selection) do s
        lc_choos[] = s
    end
    on(menupair.value) do v
        ipair[] = v
    end
    fig
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

end