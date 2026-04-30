@reexport module CubePlots
using Colors: Colorant
using Makie
using YAXArrays
using Dates, DateFormats
using ..ARCEMEAnalysis
using ..ARCEMEAnalysis: local_cubepath, arceme_legends, Event, arceme_open
using DataStructures: counter, OrderedDict
using Statistics: mean
using GeoMakie
import DimensionalData as DD
using Random: seed!

ESA_LC_colormap = [
    (color = "#006400", alpha = 255, value = 10, label = "Tree cover"),
    (color = "#ffbb22", alpha = 255, value = 20, label = "Shrubland"),
    (color = "#ffff4c", alpha = 255, value = 30, label = "Grassland"),
    (color = "#f096ff", alpha = 255, value = 40, label = "Cropland"),
    (color = "#fa0000", alpha = 255, value = 50, label = "Built-up"),
    (color = "#b4b4b4", alpha = 255, value = 60, label = "Bare / sparse vegetation"),
    (color = "#f0f0f0", alpha = 255, value = 70, label = "Snow and ice"),
    (color = "#0064c8", alpha = 255, value = 80, label = "Permanent water bodies"),
    (color = "#0096a0", alpha = 255, value = 90, label = "Herbaceous wetland"),
    (color = "#00cf75", alpha = 255, value = 95, label = "Mangroves"),
    (color = "#fae6a0", alpha = 255, value = 100, label = "Moss and lichen"),
]
clc_clrs = [
    colorant"#006400",
    colorant"#ffbb22",
    colorant"#ffff4c",
    colorant"#f096ff",
    colorant"#fa0000",
    colorant"#b4b4b4",
    colorant"#f0f0f0",
    colorant"#0064c8",
    colorant"#0096a0",
    colorant"#00cf75",
    colorant"#fae6a0",
]

cty_clrs = cgrad(:batlowS, 18, categorical = true)

ctykeys = hcat(
    values(arceme_legends["CTY"]) |> collect,
    keys(arceme_legends["CTY"]) |> collect,
    ARCEMEAnalysis.HRL.mctykeymap.(keys(arceme_legends["CTY"]) |> collect),
    ARCEMEAnalysis.HRL.ctykeymap.(keys(arceme_legends["CTY"]) |> collect),
)

function getticks(mn, mx)
    n = length(range(mn, mx))
    collect((mn+(mx-mn)/(2*n)):(mx-mn)/n:mx)
end

function plotlc(
    current_event::Event;
    batch = "ARCEME-DC-8",
    savef = true,
    fname = "fig_ESA_LC_$(current_event.uid).png",
)
    ds = arceme_open(current_event, batch = batch)
    f, ax, p = image(
        ds.ESA_LC[time = 1],
        colormap = cgrad(clc_clrs; categorical = true),
        colorrange = (10, 110),
        # label = "Land cover",
        axis = (title = "$(current_event.uid) - Land Cover", aspect = DataAspect()),
    )
    Colorbar(
        f[1, 2],
        colormap = cgrad(clc_clrs; categorical = true),
        colorrange = (1, 11),
        ticks = (getticks(1, 11), collect(values(arceme_legends["ESA_LC"]))[2:end]),
    )
    savef && save(fname, f)
    f
end

function plotcty(
    current_event::Event;
    batch = "ARCEME-DC-8",
    savef = false,
    fname = "fig_CTY_$(current_event.uid).png",
)
    ds = arceme_open(current_event, batch = batch)
    f, ax, p = image(
        ds.CTY,
        colormap = cgrad(clc_clrs; categorical = true),
        colorrange = (2, 20),
        lclip = :white,
        hclip = :white,
        axis = (title = "$(current_event.uid) - Main Crop Type", aspect = DataAspect()),
    )
    Colorbar(
        f[1, 2],
        colormap = cgrad(clc_clrs; categorical = true),
        colorrange = (2, 20),
        lclip = :white,
        hclip = :white,
        ticks = (getticks(1, 17), [k for k in ctykeys[2:18, 1]]),
    )
    savef && save(fname, f)
    f
end

yydoy2date(x::Integer) = DateTime("$(Int(2000 + x÷1e3))", "YY") + Day(x % 1e3)
yydoy2date(x::Missing) = missing
yydoy2date(x::Vector) = broadcast(yydoy2date, x)

yydoy2ydec(x::Integer) = yeardecimal(yydoy2date(x))
yydoy2ydec(x::Missing) = missing
yydoy2ydec(x::Vector) = broadcast(yydoy2ydec, x)

rep(A) = replace(
    convert(Vector{Union{Missing,Int32}}, A),
    0 => missing,
    65526 => missing,
    65527 => missing,
    65531 => missing,
    65532 => missing,
    65533 => missing,
    65534 => missing,
    65535 => missing,
)


function plotmceh(
    current_event::Event;
    batch = "ARCEME-DC-8",
    savef = true,
    fname = "fig_CPMCEH_$(current_event.uid).png",
)
    ds = arceme_open(current_event, batch = batch)
    event_date = arceme_eventdate(current_event)
    mce = ds.CPMCE.data[:] |> rep |> yydoy2ydec
    mch = ds.CPMCH.data[:] |> rep |> yydoy2ydec

    fig = Figure(size = (1000, 500))
    ax = Axis(
        fig[1, 1],
        title = current_event.uid,
        xlabel = "left: CPMCE / right: CPMCH",
        ylabel = "Year decimal",
    )
    cty = ARCEMEAnalysis.HRL.ctykeymap.(ds.CTY.data[:])[.!(ismissing.(mce))]
    clrs = cgrad(:lipariS, 14, categorical = true)
    for i = 1:14
        hist!(
            ax,
            mce[.!(ismissing.(mce))][cty.==i],
            scale_to = -0.5,
            offset = i,
            direction = :x,
            color = clrs[i],
            label = "$i: $(ctykeys[i,1])",
        )
        hist!(
            ax,
            mch[.!(ismissing.(mch))][cty.==i],
            scale_to = 0.5,
            offset = i,
            direction = :x,
            color = clrs[i],
        )
    end
    ax.xticks = 1:14
    # Colorbar(fig[1,2], label = "Main crop type",
    #     colormap = cgrad(:lipariS, 14, categorical = true),
    #     colorrange = (1,14),
    #     ticks=(getticks(1,14),ctykeys[1:14,1]),
    # )
    hlines!(
        ax,
        yeardecimal(event_date),
        color = "grey",
        label = "Event date: $(Date(event_date))",
    )
    fig[1, 2] = Legend(fig, ax, "Main Crop Types", framevisible = false)

    savef && save(fname, fig)
    fig
end

function piecty(
    current_event::Event;
    batch = "ARCEME-DC-8",
    savef = true,
    fname = "fig_CTY_pie_$(current_event.uid).png",
)
    cubename = arceme_cubename(current_event)
    fp_cty = open_dataset(
        joinpath(
            ARCEMEAnalysis.local_cubepath,
            "$(batch)-fingerprints",
            "$(split(cubename,".")[1])__fp_CTY.zarr",
        ),
    )
    f, ax, plt = pie(
        fp_cty.class_fractions.data[2:18],
        color = cty_clrs[2:end],
        radius = 4,
        inner_radius = 2,
        strokecolor = :white,
        strokewidth = 2,
        axis = (autolimitaspect = 1,),
        label = [k => (; color = c) for (k, c) in zip(ctykeys[2:18, 1], cty_clrs[2:end])],
    )
    leg = Legend(f[1, 2], ax)
    hidedecorations!(ax)

    savef && save(fname, f)
    f
end

function fpfig(
    current_event::Event,
    class;
    idx = "kNDVI",
    mclass = "Annual crop",
    batch = "ARCEME-DC-8",
)
    cubename = arceme_cubename(current_event)
    event_date = arceme_eventdate(current_event)
    fp_cty = open_dataset(
        joinpath(
            ARCEMEAnalysis.local_cubepath,
            "$(batch)-fingerprints",
            "$(split(cubename,".")[1])__fp_CTY.zarr",
        ),
    )
    mfp_cty = open_dataset(
        joinpath(
            ARCEMEAnalysis.local_cubepath,
            "$(batch)-fingerprints",
            "$(split(cubename,".")[1])__fp_MCTY.zarr",
        ),
    )
    fig = Figure()
    tempo_s2 = lookup(fp_cty.time_sentinel_2_l2a)
    tick_positions = datetime2unix.((event_date-Year(1)):Month(4):(event_date+Year(1)))  # Convert dates to Unix time (positions)
    tick_labels = [
        "$(Dates.format(t, "yyyy-mm-dd") )" for
        t = (event_date-Year(1)):Month(4):(event_date+Year(1))
    ]
    n = length(tempo_s2)
    ax = Axis(
        fig[1, 1],
        title = "$(current_event.uid)",
        xlabel = "Time",
        ylabel = idx,
        xticks = (tick_positions, tick_labels),
    )
    scatterlines!(
        ax,
        datetime2unix.(tempo_s2),
        mfp_cty.uncorrected_s2_indices[band_s2 = At(idx), class = At(mclass)].data[:],
        label = "$(mclass) - $idx",
        linewidth = 2,
    )
    scatterlines!(
        ax,
        datetime2unix.(tempo_s2),
        fp_cty.uncorrected_s2_indices[band_s2 = At(idx), class = At(class)].data[:],
        label = "$(class) - $idx",
        linewidth = 1,
    )
    # ax.xticks = ([1,n÷2+1,n], string.(Date.(tempo_s2[[1,n÷2+1,n]])))
    fig[2, 1] = Legend(fig, ax, framevisible = false, orientation = :horizontal, nbanks = 1)
    fig
end

function getsplines(current_event; bands = ["kNDVI", "NDVI"], batch = "ARCEME-DC-8")
    ds = arceme_open(current_event, batch = batch)
    firstday = arceme_eventdate(current_event) - Year(1)
    lastday = arceme_eventdate(current_event) + Year(1)
    nms = Millisecond(lastday - firstday).value
    linscal = d -> Millisecond(d - firstday).value * 2 / nms
    x_s2 = map(linscal, ds.time_sentinel_2_l2a.val)
    outax = YAXArrays.time(range(firstday, lastday, step = Day(1)))
    s2 = YAXArray(
        (ds.x, ds.y, ds.time_sentinel_2_l2a, Dim{:band_s2}(bands)),
        # !!cat(([ds[bands[i]].data] for i in eachindex(bands)), dims = 4),
        cat(ds.kNDVI.data, ds.NDVI.data, dims = 4),
    )
    r_s2 = xmap(
        ARCEMEAnalysis.FingerprintPlots.filter_spline,
        s2 ⊘ :time_sentinel_2_l2a,
        output = XOutput(outax),
        function_args = (x_s2, length(outax)),
    )
    return r_s2
end

function plot_cpbsb_splines(
    current_event::Event,
    crop;
    batch = "ARCEME-DC-8",
    idx = "kNDVI",
    nsamples = 100,
    nmonths = 6,
    random = true,
    ind = [1, 10, 100, 1000],
)
    ds = arceme_open(current_event, batch = batch)
    r_s2 = getsplines(current_event)
    event_date = arceme_eventdate(current_event)

    cls = ctykeys[findfirst(ctykeys[:, 1] .== crop), 2]
    indcrop = findall(x -> x == cls, ds.CTY.data[:, :])
    f = Figure()
    tick_positions =
        datetime2unix.((event_date-Year(1)):Month(nmonths):(event_date+Year(1)))  # Convert dates to Unix time (positions)
    tick_labels = [
        "$(Dates.format(t, "yyyy-mm-dd") )" for
        t = (event_date-Year(1)):Month(nmonths):(event_date+Year(1))
    ]
    a1 = Axis(
        f[1, 1],
        title = "$(current_event.uid) - $crop",
        xlabel = "Date",
        ylabel = idx,
        xticks = (tick_positions, tick_labels),
    )
    if random
        Random.seed!(42)
        cind = rand(indcrop, nsamples)
    else
        cind = indcrop[ind]
    end
    for i in cind
        cpbsb = ds.CPBSB.data[i[1], i[2]]
        cpbsb > 65000 && continue
        lines!(
            a1,
            datetime2unix.(lookup(r_s2.time)),
            r_s2[band_s2 = At(idx)].data[:, i[1], i[2], 1],
            color = cpbsb,
            colorrange = (0, 295),
            colormap = (:viridis, 0.5),
        )
    end
    Colorbar(
        f[1, 2],
        colorrange = (0, 295),
        colormap = (:viridis, 0.5),
        label = "Bare soil before crop emergence (days)",
    )
    f
end


function crop_splines(
    current_event,
    crop;
    mx = 0.25,
    nsamples = 100,
    batch = "ARCEME-DC-8",
    idx = "kNDVI",
)
    ds = arceme_open(current_event, batch = batch)
    event_date = arceme_eventdate(current_event)
    r_s2 = getsplines(current_event, batch = batch)

    cls = ctykeys[findfirst(ctykeys[:, 1] .== crop), 2]
    indcrop = findall(x -> x == cls, ds.CTY.data[:, :])

    f = Figure()
    tempo = datetime2unix.(lookup(r_s2.time))
    tick_positions = datetime2unix.((event_date-Year(1)):Month(6):(event_date+Year(1)))  # Convert dates to Unix time (positions)
    tick_labels = [
        "$(Dates.format(t, "yyyy-mm-dd") )" for
        t = (event_date-Year(1)):Month(6):(event_date+Year(1))
    ]
    a1 = Axis(
        f[1, 1],
        title = "$(current_event.uid) - $(crop)",
        ylabel = idx,
        xticks = (tick_positions, tick_labels),
    )
    Random.seed!(42)
    for i in rand(indcrop, nsamples)
        ts = ds[idx].data[i[1], i[2], :]
        mn = minimum(ts[.!isnan.(ts)])
        lines!(
            a1,
            tempo,
            r_s2[band_s2 = At(idx)].data[:, i[1], i[2], 1],
            color = mn,
            colorrange = (0, mx),
            colormap = (:viridis, mx),
        )
    end
    Colorbar(
        f[1, 2],
        colorrange = (0, mx),
        colormap = (:viridis, mx),
        label = "Minimum $idx",
    )
    f
end


export plotlc, plotmceh, piecty, fpfig, getsplines, plot_cpbsb_splines
end