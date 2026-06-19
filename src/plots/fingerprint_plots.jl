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
using DiskArrays: cache
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



function most_common_class(pairnow)
    ds1, ds2 = arceme_open.(pairnow)
    #Find maximum landcover in both cubes
    classfrac1, classfrac2 = (ds1.class_fractions[:], ds2.class_fractions[:])
    (_, i1), (_, i2) = findmax.((classfrac1, classfrac2))
    if classfrac1[i1] * classfrac2[i1] > classfrac1[i2] * classfrac2[i2]
        ds1.class.val[i1]
    else
        ds1.class.val[i2]
    end
end

function pairwise_without_missings(x, y)
    out = Point2f[]
    indices = Int[]
    for i in eachindex(x, y)
        if !ismissing(y[i]) && isfinite(y[i])
            push!(out, Point2f(x[i], y[i]))
            push!(indices, i)
        end
    end
    out, indices
end

function open_plot_data(ds, class; strata="ESA_LC")
    firstday = arceme_eventdate(ds) - Year(1)
    lastday = arceme_eventdate(ds) + Year(1)
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
    index_orig = OrderedDict()

    for band in ds.band_s2
        interplots[band], _ = pairwise_without_missings(range(0, 2, length=length(outax)), r_s2[band_s2=DD.At(band)].data[:])
        dataplots[band], index_orig[band] = pairwise_without_missings(x_s2, s2[band_s2=DD.At(band)].data[:])
    end

    for band in ds.band_s1
        interplots[band], _ = pairwise_without_missings(range(0, 2, length=length(outax)), r_s1[band_s1=DD.At(band)].data[:])
        dataplots[band], index_orig[band] = pairwise_without_missings(x_s1, s1[band_s1=DD.At(band)].data[:])
    end

    classfracs = ds.class_fractions.data[:]
    sfracs = sortperm(classfracs, rev=true)
    strs = map(sfracs[1:3]) do ifrac
        string(collect(values(arceme_legends[strata]))[ifrac], " => ", round(classfracs[ifrac] * 100), "%")
    end
    classdesc = join(strs, ", ")


    (dataplots, interplots, classdesc, index_orig)
end


function fingerprint_plot(pairid=1; class=nothing, plotdict=IdDict(), interactive=true, strata="ESA_LC", mask_clouds=false)

    ipair = Observable{Int}(pairid)
    alleventpairs = arceme_validpairs()

    if class === nothing
        class = most_common_class(alleventpairs[pairid])
    end

    class_choos = Observable{Any}(class)

    ds = lift(ipair) do i
        ev = alleventpairs[i]
        cache.(arceme_open.(ev))
    end

    d = lift(ds, class_choos) do dss, class
        open_plot_data(dss[1], class)
    end
    dhp = lift(ds, class_choos) do dss, class
        open_plot_data(dss[2], class)
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

    repr_ind_d = Observable{Any}(nothing)
    repr_ind_dhp = Observable{Any}(nothing)
    omask_clouds = Observable{Bool}(mask_clouds)

    npl = maximum(first, keys(allax))
    for ipl in 1:npl
        linkyaxes!(allax[(ipl, 1)], allax[(ipl, 2)])
    end
    ov1 = Axis(fig[5:6, 3], aspect=DataAspect(), title="Drought")
    img1 = lift((y, i, ind, mask_clouds) -> arceme_representative_image(i[1]; class=y, repr_ind=ind, mask_clouds), class_choos, ds, repr_ind_d, omask_clouds)
    heatmap!(ov1, img1)
    ov2 = Axis(fig[7:8, 3], aspect=DataAspect(), title="Drought + HP")
    img2 = lift((y, i, ind, mask_clouds) -> arceme_representative_image(i[2]; class=y, repr_ind=ind, mask_clouds), class_choos, ds, repr_ind_dhp, omask_clouds)
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

        #menupair = Slider(fig, range=1:length(alleventpairs), startvalue=pairid, update_while_dragging=false)

        menupair = Textbox(fig, placeholder=string(ipair[]), validator=Int)
        menuclass = Menu(fig, options=collect(values(arceme_legends[strata])), default=class)
        # cb = Checkbox(checkboxes[1, 1], checked=false)
        cb2 = Checkbox(fig, checked=true)
        # cb3 = Checkbox(checkboxes[3, 1], checked=true)
        on(c->setindex!(omask_clouds, c), cb2.checked)

        fig[1:2, 3] = vgrid!(
            menupair,
            Label(fig, lift(i -> i[3], d), color=:red),
            Label(fig, lift(i -> i[3], dhp), color=:blue),
            menuclass,
            cb2,
        )


        #fig[1:2,3] = allmenus


        # checkboxes = GridLayout(allmenus[5:7, 1:3])




        # Label(checkboxes[1, 2], "Classmask", halign=:left)
        # Label(checkboxes[2, 2], "Cloudmask", halign=:left)
        # Label(checkboxes[3, 2], "SCL Mask", halign=:left)
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
        on(events(fig).mousebutton) do event
            if event.button == Mouse.left && event.action == Mouse.press
                p, i = pick(fig)
                index, i_dhp = get(plotdict, p, (nothing, nothing))
                if !isnothing(index)
                    if i_dhp == 1
                        indices = last(d[])
                        repr_ind_d[] = indices[index][i]
                    else
                        indices = last(dhp[])
                        repr_ind_dhp[] = indices[index][i]
                    end
                    Consume(true)
                end
            end
        end
    end
    fig
end


arceme_representative_image(ev::Event; index_use="NDVI", brighten_factor=2, repr_ind=nothing, class=nothing, mask_clouds=false) =
    arceme_representative_image(arceme_open(ev); index_use, brighten_factor, repr_ind, class, mask_clouds)

"""
    arceme_representative_image(ev)

Creates an RGB image highlighting on a cloud-free time step with high NDVI for a certain land cover type 
"""
function arceme_representative_image(ds; index_use="NDVI", brighten_factor=2, repr_ind=nothing, class=nothing, strata="ESA_LC", mask_clouds=false)
    if class === nothing
        classfrac = ds.class_fractions.data[:]
        _, iclass = findmax(classfrac)
        class = ds.class.val[iclass]
    end
    if repr_ind === nothing
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
    else
        @show repr_ind
        @show ds.time_sentinel_2_l2a[repr_ind]
    end
    brfilt(a, factor, alpha) = RGBA(clamp(a.r * factor, 0, 1), clamp(a.g * factor, 0, 1), clamp(a.b * factor, 0, 1), alpha)
    ilc = findfirst(==(class), collect(values(arceme_legends[strata])))
    r = broadcast(arceme_rgb(ds)[time_sentinel_2_l2a=repr_ind], ds[strata][:, :, 1], brighten_factor) do col, lc, f
        lcv = ARCEMEAnalysis.lckeymap(lc)
        alpha = lcv == ilc ? 1.0 : 0.5
        brfilt(col, f, alpha)
    end
    if mask_clouds
        cl = ds.cloud_mask[:, :, repr_ind]
        scl = ds.SCL[:, :, repr_ind]
        is_cloud = map(ARCEMEAnalysis._is_cloud, cl.data, scl.data)
        res = r[:, :]
        res[is_cloud] .= RGBA(0.0, 0.0, 1.0, 1.0)
        res
    else
        r[:, :]
    end
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