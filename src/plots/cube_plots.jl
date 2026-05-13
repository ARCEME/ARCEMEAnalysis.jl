@reexport module CubePlots
using Colors: Colorant
using Makie
using YAXArrays
using Dates, DateFormats
using ..ARCEMEAnalysis
using ..ARCEMEAnalysis: local_cubepath, arceme_legends, Event, arceme_open
using DataStructures: counter, OrderedDict
using Statistics: mean, median
using GeoMakie
import DimensionalData as DD
using Random: seed!
using DataFrames

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

cty_clrs = cgrad(
    vcat(
        cgrad(:lipariS, 14, categorical = true).colors.colors,
        cgrad(:bamakoS, 4, categorical = true).colors.colors,
    ),
    18,
    categorical = true,
)

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
    hidedecorations!(ax)
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
    savef = true,
    fname = "fig_CTY_$(current_event.uid).png",
)
    ds = arceme_open(current_event, batch = batch)
    cty = map(ARCEMEAnalysis.HRL.ctykeymap, ds.CTY)
    f, ax, p = image(
        cty,
        colormap = cgrad(cty_clrs; categorical = true),
        colorrange = (1, 18),
        # lowclip = :white,
        highclip = :white,
        axis = (title = "$(current_event.uid) - Main Crop Type", aspect = DataAspect()),
    )
    hidedecorations!(ax)
    Colorbar(
        f[1, 2],
        colormap = cgrad(cty_clrs; categorical = true),
        colorrange = (1, 18),
        # lowclip = :white,
        highclip = :white,
        ticks = (getticks(1, 18), [k for k in ctykeys[1:18, 1]]),
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
    for i = 1:14
        hist!(
            ax,
            mce[.!(ismissing.(mce))][cty.==i],
            scale_to = -0.5,
            offset = i,
            direction = :x,
            color = cty_clrs.colors.colors[i],
            label = "$i: $(ctykeys[i,1])",
        )
        hist!(
            ax,
            mch[.!(ismissing.(mch))][cty.==i],
            scale_to = 0.5,
            offset = i,
            direction = :x,
            color = cty_clrs.colors.colors[i],
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
    savef = true,
    fname = "fig_fp_$(class)_$(current_event.uid).png",
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
    if !isnothing(mclass)
        mfp_cty = open_dataset(
            joinpath(
                ARCEMEAnalysis.local_cubepath,
                "$(batch)-fingerprints",
                "$(split(cubename,".")[1])__fp_MCTY.zarr",
            ),
        )
        scatterlines!(
            ax,
            datetime2unix.(tempo_s2),
            mfp_cty.uncorrected_s2_indices[band_s2 = At(idx), class = At(mclass)].data[:],
            label = "$(mclass) - $idx",
            linewidth = 2,
        )
    end
    scatterlines!(
        ax,
        datetime2unix.(tempo_s2),
        fp_cty.uncorrected_s2_indices[band_s2 = At(idx), class = At(class)].data[:],
        label = "$(class) - $idx",
        linewidth = 1,
    )
    # ax.xticks = ([1,n÷2+1,n], string.(Date.(tempo_s2[[1,n÷2+1,n]])))
    fig[2, 1] = Legend(fig, ax, framevisible = false, orientation = :horizontal, nbanks = 1)

    savef && save(fname, fig)
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
    savef = true,
    fname = "fig_$(idx)_splines_$(crop)_cpbsb_$(current_event.uid).png",
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
        seed!(42)
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
    savef && save(fname, f)
    f
end

function getindcrop(current_event, crop; batch = "ARCEME-DC-8")
    ds = arceme_open(current_event, batch = batch)
    cls = ctykeys[findfirst(ctykeys[:, 1] .== crop), 2]
    indcrop = findall(x -> x == cls, ds.CTY.data[:, :])
end

function crop_splines(
    current_event,
    crop;
    mx = 0.25,
    nsamples = 100,
    batch = "ARCEME-DC-8",
    idx = "kNDVI",
)
    ds = arceme_open(current_event; batch)
    event_date = arceme_eventdate(current_event)
    r_s2 = getsplines(current_event; batch)

    indcrop = getindcrop(current_event, crop; batch)

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
    seed!(42)
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

function tsstatsA(xout, ts, mce, mch, tempo)
    if (ismissing(mce) || ismissing(mch))
        xout .= missing
        return nothing
    end
    missmask = findall(x -> !ismissing(x) && isfinite(x), ts)
    ind = tempo[missmask] .> mce .&& tempo[missmask] .< mch
    mx = maximum(ts[missmask][ind])
    # tmx = Day(tempo[missmask][ind][ts[missmask][ind].==mx][1] - tempo[1]).value
    tmx = dayofyear(tempo[missmask][ind][ts[missmask][ind].==mx][1])
    sts = sum(ts[missmask][ind])
    xout .= cat(mx, tmx, sts, dims = 1)
    return nothing
end

"""
    
"""
function arceme_acrop(
    current_event,
    crop;
    batch = "ARCEME-DC-8",
    savef = true,
    fname = "fig_violin_anova_$(crop)_$(current_event.uid).png",
    printres = true,
)
    ds = arceme_open(current_event; batch)

    indcrop = getindcrop(current_event, crop; batch)

    gbsb_labels = ["<35", "35-60", ">60"]
    gbsb = map(ds.CPBSB) do x
        (x == 0 || x >= 65500) && return missing
        x < 35 && return 1
        (x >= 35 || x < 60) && return 2
        x >= 60 && return 3
    end

    mce = map(ds.CPMCE) do x
        (x == 0 || x >= 65500) && return missing
        return CubePlots.yydoy2date(x)
    end

    mch = map(ds.CPMCH) do x
        (x == 0 || x >= 65500) && return missing
        return CubePlots.yydoy2date(x)
    end

    outax = Dim{:mcstats}(["mx", "tmx", "sum"])
    mcstata = xmap(
        tsstatsA,
        ds.kNDVI ⊘ :time_sentinel_2_l2a,
        mce,
        mch,
        output = XOutput(outax, outtype = Union{Missing,Float32}),
        function_args = (lookup(ds.time_sentinel_2_l2a),),
    )
    data_anova = hcat(
        permutedims(mcstata.data[:, indcrop[indcrop.∈(indBSB,)], 1]),
        gbsb.data[indcrop[indcrop.∈(indBSB,)]],
    )
    df_anova =
        DataFrame(
            mx = data_anova[:, 1],
            tmx = data_anova[:, 2],
            sumkndvi = data_anova[:, 3],
            gbsb = data_anova[:, 4],
        ) |> (df -> dropmissing!(df))
    gdf = groupby(df_anova, :gbsb)
    means = gdf |> (gdf -> combine(gdf, :mx => mean, :tmx => mean, :sumkndvi => mean))
    ylabs = ["max(kNDVI)", "DOY | max(kNDVI)", "∑(kNDVI)"]

    f = Figure(size = (400, 800))
    for i = 1:3
        ax = Axis(
            f[i, 1],
            xlabel = "Bare Soil Before Main Crop emergence (Days)",
            xticks = (1:nrow(means), gbsb_labels[1:nrow(means)]),
            ylabel = ylabs[i],
        )
        violin!(ax, df_anova[:, 4], df_anova[:, i], scale = :count, show_median = true)
        for j = 1:nrow(means)
            hlines!(
                ax,
                means[j, i+1],
                xmin = (j - 1) / nrow(means),
                xmax = j / (nrow(means)),
                color = :grey,
            )
        end
        i < 3 && hidexdecorations!(ax, grid = false)
    end
    Legend(
        f[4, 1],
        [
            LineElement(color = :black, linestyle = nothing),
            LineElement(color = :grey, linestyle = nothing),
        ],
        ["Median", "Mean"],
        orientation = :horizontal,
    )
    f
    savef && save(fname, f)
    res = Dict()
    for istat in names(df_anova)[1:3]
        printres && println(istat)
        res[istat] = OneWayANOVATest(gdf[1][!, Symbol(istat)], gdf[2][!, Symbol(istat)])
        printres && println(res[istat])
    end
    return (f, res, data_anova)
end

function tsstatsP(xout, ts, tempo)
    event_date = tempo[length(tempo)÷2+1]
    missmask = findall(x -> !ismissing(x) && isfinite(x), ts)
    indWinter = month.(tempo[missmask]) .∈ ([1, 2, 11, 12],)
    indGrowing1 = month.(tempo[missmask]) .∈ (3:10,) .&& tempo[missmask] .< event_date
    indGrowing2 = month.(tempo[missmask]) .∈ (3:10,) .&& tempo[missmask] .>= event_date
    avr = mean(ts[missmask][indWinter])
    mx1 = maximum(ts[missmask][indGrowing1])
    tmx1 = dayofyear(tempo[missmask][indGrowing1][ts[missmask][indGrowing1].==mx1][1])
    sts1 = sum(ts[missmask][indGrowing1])
    mx2 = maximum(ts[missmask][indGrowing2])
    tmx2 = dayofyear(tempo[missmask][indGrowing2][ts[missmask][indGrowing2].==mx2][1])
    sts2 = sum(ts[missmask][indGrowing2])
    xout .= cat(avr, mx1, tmx1, sts1, mx2, tmx2, sts2, dims = 1)
    return nothing
end


function groupwinterkndvi(x, med)
    x < med && return 1
    x >= med && return 2
end

function arceme_pcrop(
    current_event,
    crop;
    batch = "ARCEME-DC-8",
    savef = true,
    fname = "fig_violin_anova_$(crop)_$(current_event.uid).png",
    printres = true,
)
    ds = arceme_open(current_event; batch)
    outax = Dim{:mcstats}(["avr", "mx1", "tmx1", "sum1", "mx2", "tmx2", "sum2"])
    grapestata = xmap(
        tsstatsP,
        ds.kNDVI ⊘ :time_sentinel_2_l2a,
        output = XOutput(outax, outtype = Union{Missing,Float32}),
        function_args = (lookup(ds.time_sentinel_2_l2a),),
    )
    indcrop = CubePlots.getindcrop(current_event, crop, batch = "ARCEME-DC-6")
    data_anova = permutedims(grapestata.data[:, indcrop, 1])
    df_anova =
        DataFrame(
            avr = data_anova[:, 1],
            mx1 = data_anova[:, 2],
            tmx1 = data_anova[:, 3],
            sumkndvi1 = data_anova[:, 4],
            mx2 = data_anova[:, 5],
            tmx2 = data_anova[:, 6],
            sumkndvi2 = data_anova[:, 7],
        ) |> (df -> dropmissing!(df))
    med = median(df_anova.avr) # => balanced samples
    gwinter_labels = ["< median", ">= median"]
    transform!(df_anova, :avr => ByRow(x -> groupwinterkndvi(x, med)) => :gwinter)
    gdf = groupby(df_anova, :gwinter)
    # compute mean over DJF
    means =
        gdf |> (
            gdf -> combine(
                gdf,
                :mx1 => mean,
                :tmx1 => mean,
                :sumkndvi1 => mean,
                :mx2 => mean,
                :tmx2 => mean,
                :sumkndvi2 => mean,
            )
        )
    ylabs = [
        "max(kNDVI) | Year 1",
        "DOY | max(kNDVI) | Year 1",
        "∑(kNDVI) | Year 1",
        "max(kNDVI) | Year 2",
        "DOY | max(kNDVI) | Year 2",
        "∑(kNDVI) | Year 2",
    ]
    f = Figure(size = (600, 800))
    for i = 1:6
        ax = Axis(
            f[i > 3 ? i - 3 : i, i > 3 ? 2 : 1],
            xlabel = "Average winter kNDVI",
            xticks = (1:nrow(means), gwinter_labels[1:nrow(means)]),
            ylabel = ylabs[i],
        )
        violin!(ax, df_anova[:, 8], df_anova[:, i+1], scale = :count, show_median = true)
        for j = 1:nrow(means)
            hlines!(
                ax,
                means[j, i+1],
                xmin = (j - 1) / nrow(means),
                xmax = j / (nrow(means)),
                color = :grey,
            )
        end
        !(i ∈ [3, 6]) && hidexdecorations!(ax, grid = false)
    end
    Legend(
        f[4, 1:2],
        [
            LineElement(color = :black, linestyle = nothing),
            LineElement(color = :grey, linestyle = nothing),
        ],
        ["Median", "Mean"],
        orientation = :horizontal,
    )
    savef && save(fname, f)
    # ANOVA
    res = Dict()
    for istat in names(df_anova)[2:7]
        printres && println(istat)
        res[istat] = OneWayANOVATest(gdf[1][!, Symbol(istat)], gdf[2][!, Symbol(istat)])
        printres && println(res[istat])
    end
    return (f, res, data_anova)
end

"""
`plot_ts(idx::CartesianIndex, current_event::Event; layer="kNDVI")`

plot time series.
return (f, tmp1, tmp_rm, sp)
"""
function plot_ts(idx, current_event; layer = "kNDVI", batch = "ARCEME-DC-8")
    ds = arceme_open(current_event, batch = batch)
    tempo = lookup(ds.time_sentinel_2_l2a)
    dayssincestart = datetime2julian.(tempo) .- datetime2julian(tempo[1])
    tmp = ds[layer].data[idx[1], idx[2], :]
    tick_positions =
        datetime2julian.((event_date-Year(1)):Month(6):(event_date+Year(1))) .-
        datetime2julian(tempo[1])  # Convert dates to Unix time (positions)
    tick_labels = [
        "$(Dates.format(t, "yyyy-mm-dd") )" for
        t = (event_date-Year(1)):Month(6):(event_date+Year(1))
    ]
    f = Figure()
    a = Axis(f[1, 1], xticks = (tick_positions, tick_labels))
    p = scatter!(a, dayssincestart, tmp, label = "Cloud filtered")
    missmask = findall(x -> !ismissing(x) && isfinite(x), tmp)
    # flag jumps
    indflag = findall(
        broadcast(
            (x, y) -> x > 0.5 && y <= 20,
            diff(tmp[missmask]),
            diff(dayssincestart[missmask]),
        ),
    )
    deleteat!(missmask, indflag .+ 1)
    # scatterlines!(a, dayssincestart[missmask], tmp[missmask], label = "Jump filtered")
    # linear interpolation 
    tmp1 = fill(NaN, length(tmp))
    tmp1[missmask] .= tmp[missmask]
    nodes = (dayssincestart[missmask],)
    itp = interpolate(nodes, tmp[missmask], Gridded(Interpolations.Linear()))
    if missmask[1] > 1
        tmp1[1:(missmask[1]-1)] .= tmp[missmask[1]]
    end
    if missmask[end] < length(dayssincestart)
        tmp1[(missmask[end]+1):end] .= tmp[missmask[end]]
    end
    tmp1[missmask[1]:missmask[end]] = itp(dayssincestart[missmask[1]:missmask[end]])

    scatterlines!(
        a,
        dayssincestart,
        tmp1,
        color = Makie.wong_colors()[1],
        markersize = 5,
        label = "Linear interpolation",
    )
    f
    # # 2. SavitzkyGolay filter 
    # tmp_sg1 = savitzky_golay(tmp1, 7, 2)
    # lines!(a, dayssincestart, tmp_sg1.y, label = "Savitzky–Golay w=3, o=1")
    # # tmp_sg2 = savitzky_golay(tmp_sg1.y, 5, 3)
    # # lines!(dayssincestart, tmp_sg2.y, label = "Savitzky–Golay w=5, o=3")
    # # tmp_sg3 = savitzky_golay(tmp_sg2.y, 7, 2)
    # # lines!(dayssincestart, tmp_sg3.y, label = "Savitzky–Golay w=7, o=5")
    # # tmp_sg4 = savitzky_golay(tmp_sg3.y, 9, 6)
    # # lines!(dayssincestart, tmp_sg4.y, label = "Savitzky–Golay w=9, o=6")
    # f
    # # the succesive iterations don't do anything! Because of the linear interpolation, I suppose?
    # # Maybe I should use a rolling mean before

    # rolling mean
    tmp_rm = running(mean, tmp1, 5)
    lines!(
        a,
        dayssincestart[1:end-2],
        tmp_rm[3:end],
        color = Makie.wong_colors()[2],
        label = "Rolling mean w=5",
    )
    r_s2 = CubePlots.getsplines(current_event)
    sp = r_s2[band_s2 = At(layer)].data[:, idx[1], idx[2], 1]
    lines!(a, sp, color = Makie.wong_colors()[3], label = "Spline 19 knots")
    axislegend(a)
    f
    return (f, tmp1, tmp_rm, sp)
end

export plotlc,
    plotcty,
    plotmceh,
    piecty,
    fpfig,
    getsplines,
    plot_cpbsb_splines,
    crop_splines,
    getindcrop,
    arceme_pcrop,
    plot_ts
end