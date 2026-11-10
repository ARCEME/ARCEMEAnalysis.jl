# recovery.jl
#
# Recovery analysis: compare vegetation recovery after drought between
# drought-only (d) and drought + heavy precipitation (dhp) event pairs, by
# averaging aligned spline-fitted index time series across data-cube pairs.
#
# The analysis lives in the `RecoveryAnalysis` submodule so its (many) helpers
# do not pollute the `ARCEMEAnalysis` namespace. Public entry points are the
# `arceme_`-prefixed functions listed in the export statement; everything else
# is private to the submodule.

module RecoveryAnalysis

import ..ARCEMEAnalysis
using ..ARCEMEAnalysis
using ..ARCEMEAnalysis.FingerprintPlots: filter_spline
using YAXArrays
using Dates
using Statistics: quantile, mean, std
using StatsBase: sem
import DimensionalData as DD
import HypothesisTests
using Makie
import CSV

export arceme_recovery_cache, arceme_recovery_splines, arceme_recovery_pairs,
    arceme_recovery_test, arceme_recovery_plot, arceme_recovery_plot_delta,
    arceme_recovery_export_csv, arceme_run_recovery,
    VegStateGrouper, GrowingSeasonSplitter, LandCoverSelector,
    vegstate_group, vegstate_groups, group_token

# Event date lands at this index on the 2-year daily time_smooth axis
# (eventdate ± 365 days => 731 daily points, day 0 = eventdate at index 366),
# fixed for every cube so curves align across pairs by index.
const EVENT_INDEX = 366

# Shared time axis used ONLY when stacking many events, expressed as days since
# the event onset (index 366 = 0). compute_indexsplines itself stays in real
# date space; this relative axis is applied at aggregation time via setdims to
# (a) destroy per-pair date info so every site's axis is identical (required by
# concatenatecubes, which errors on differing axes) and (b) label the axis with
# the relative-to-event semantics the veg-state grouper and recovery window
# stratify on. Length 731, same order as the date grid, so index 366 = day 0.
const REL_TIME_AXIS = DD.Dim{:time_smooth}(-365:365)

# Minimum number of sites with data before a cell's cross-site reduction is
# reported; below this the cell is left `missing`.
const MIN_SAMPLES = 3

# FDR significance level for the per-day paired tests
const FDR_ALPHA = 0.05

const S2_BANDS = ("NDVI", "kNDVI", "NDWI", "EVI2", "NIRv", "NDMI", "NSDSI3", "WDRVI")
const S1_BANDS = ("DpRVIVV", "vv_db", "vh_db")

# ---------------------------------------------------------------------------
# Vegetation-state groupers
# ---------------------------------------------------------------------------

abstract type VegStateGrouper end

"""
    vegstate_group(g::VegStateGrouper, curve_d, curve_h, i_ev)

Classify the selected land-cover class into a vegetation-state group for a cube
pair, reconciling the d and dhp cubes. `curve_d`/`curve_h` are the grouper-index
curves (length-731 vectors, one per cube) for the selected LC; `i_ev` is the
event-date index on the time axis.

Returns a `Symbol` from the grouper's group set (see `vegstate_groups`), or
`:skip` on genuine d/dhp disagreement.
"""
function vegstate_group end

"""
    vegstate_groups(g::VegStateGrouper) -> Tuple{Symbol,...}

The set of veg-state group `Symbol`s this grouper can emit (excluding `:skip`).
`arceme_run_recovery` iterates this to build the land-cover × group loop.
Defaults to `(:growing_season, :dormant)`; a grouper overrides it to declare
its own groups.
"""
vegstate_groups(::VegStateGrouper) = (:growing_season, :dormant)

"""
    group_token(g::VegStateGrouper, group::Symbol) -> String

Output-directory token (folder name) for a veg-state group `Symbol` under this
grouper's output tree. Defaults to `string(group)`; a grouper overrides it for
human-readable folder names.
"""
group_token(g::VegStateGrouper, group::Symbol) = string(group)

"""
    GrowingSeasonSplitter <: VegStateGrouper

Threshold classifier: a class is `:growing_season` at the event when its grouper
index spline value there exceeds `vlow = quantile(index_spline, 0.95) / 3` (q95
taken over the full 2-year window). The grouper index is kNDVI when present,
falling back to NDVI for cubes lacking the kNDVI band. Reconciliation follows
RECOVERY_PLAN.md decision 1: agree -> accept; disagree but a cube is within
`tol*vlow` of its threshold -> take the d cube's label; otherwise -> `:skip`.
"""
struct GrowingSeasonSplitter <: VegStateGrouper
    tol::Float64
    threshold::Float64
end
GrowingSeasonSplitter(; tol=0.05, threshold=1/3) = GrowingSeasonSplitter(tol, threshold)

vegstate_groups(::GrowingSeasonSplitter) = (:growing_season, :dormant)
function group_token(::GrowingSeasonSplitter, group::Symbol)
    group === :growing_season && return "growing"
    group === :dormant && return "dormant"
    throw(ArgumentError("no output-directory token for veg-state group $group"))
end

Base.string(g::GrowingSeasonSplitter) = "GrowingSeasonSplitter_$(g.threshold)"

function Base.string(g::VegStateGrouper)
    string(nameof(typeof(g)))
end

# Filename token for the sensor a recovery file describes. Only Sentinel-2 is
# aggregated today; Sentinel-1 will add `-s1` sibling files without changing
# the layout.
_sensor_token() = "s2"

# Values that are neither missing nor NaN (the spline can emit NaN on fit
# failure, which skipmissing does not drop).
_valid_vals(curve) = (v for v in curve if !ismissing(v) && !isnan(v))

# q95 of a curve, skipping missing AND NaN; NaN when the curve has no valid
# values (caller treats it as "no data for this class" -> :skip).
function _gq95(curve)
    vals = collect(_valid_vals(curve))
    isempty(vals) && return NaN
    quantile(vals, 0.95)
end

# Event-point value, or `nothing` when missing.
_eventval(curve, i_ev) = ismissing(curve[i_ev]) ? nothing : curve[i_ev]

function _reconcile_growing(kd, kh, i_ev; tol, threshold)
    qd, qh = _gq95(kd), _gq95(kh)
    (isnan(qd) || isnan(qh)) && return :skip   # a cube has no valid data for this class
    vd, vh = qd * threshold, qh * threshold
    ed, eh = _eventval(kd, i_ev), _eventval(kh, i_ev)
    (ed === nothing || eh === nothing) && return :skip
    gd = ed > vd
    gh = eh > vh
    if gd == gh
        return gd ? :growing_season : :dormant
    end
    near_d = abs(ed - vd) < tol * vd
    near_h = abs(eh - vh) < tol * vh
    if near_d || near_h
        return gd ? :growing_season : :dormant
    end
    :skip
end

# Veg-state for the selected LC: run the reconciliation on the two grouper
# curves (one per cube). The curves carry `missing` (and possibly `NaN`) where
# no data; the grouper reads q95 via `_valid_vals` and skips when either cube
# has no valid values for that LC (see _reconcile_growing).
function vegstate_group(g::GrowingSeasonSplitter, curve_d, curve_h, i_ev)
    (isempty(_valid_vals(curve_d)) || isempty(_valid_vals(curve_h))) && return :skip
    _reconcile_growing(curve_d, curve_h, i_ev; tol=g.tol, threshold=g.threshold)
end

# ---------------------------------------------------------------------------
# Land-cover selection (thin, applied AFTER splines)
# ---------------------------------------------------------------------------

"""
    LandCoverSelector

Selects which land-cover class(es) to aggregate for one pass of the outer LC
loop. For now it holds a single `class` label string; the future merge step
(e.g. Grassland OR Shrubland, fractions-weighted) slots in as another method
without touching the aggregation.
"""
struct LandCoverSelector
    class::String
end
LandCoverSelector(; class) = LandCoverSelector(class)

# True when the class is present in the cube's class axis.
present(ds, s::LandCoverSelector) = s.class in ds.class.val

# Grouper-index curve for the selected class from an S2 spline cube, for the
# given grouper index `gidx` ("kNDVI" or "NDVI"). The grouper index is chosen
# ONCE per pair (symmetric across d and dhp) so the reconciliation never
# compares different indices. Returns a length-731 vector with `missing` intact.
_grouper_curve(s2, cls, gidx) = select_lc(s2, cls)[band_s2 = DD.At(gidx)].data

"""
    select_lc(s2, cls)

Select the land-cover class `cls` from an (expanded) spline cube, returning the
`[time_smooth × band_s2]` (or `[time_smooth × band_s1]`) YAXArray slice. Pure
slice; no Matrix/NaN conversion. Caller is responsible for checking the class is
present first.
"""
select_lc(s2, cls) = s2[class=DD.At(cls)]

# Selected class fraction (metadata only, not used in aggregation).
function class_fraction(ds, cls)
    icls = findfirst(==(cls), ds.class.val)
    icls === nothing ? missing : Float64(ds.class_fractions.data[icls])
end

# ---------------------------------------------------------------------------
# Spline fitting (dimension-agnostic, adapted from global_analysis.jl)
# ---------------------------------------------------------------------------

function compute_indexsplines(ds)
    # Fixed 731-day grid centered on the event date (day 0 = event date at
    # index 366) for every cube. Using `± Year(1)` would yield 731 or 732
    # points depending on leap-year crossings, misaligning day 0 across pairs;
    # `± Day(365)` is always exactly 731 and identical for all cubes.
    firstday = ARCEMEAnalysis.arceme_eventdate(ds) - Day(365)
    lastday = ARCEMEAnalysis.arceme_eventdate(ds) + Day(365)
    outax = DD.format(DD.Dim{:time_smooth}(range(firstday, lastday, step = Day(1))))
    nms = Millisecond(lastday - firstday).value
    linscal = d -> Millisecond(d - firstday).value * 2 / nms

    # xmap broadcasts filter_spline over every axis other than the time axis
    # (including the class axis), so this is dimension-agnostic: it only needs
    # the S1/S2 time axes to exist. LC selection happens downstream in select_lc.
    s1splines = if :s1_indices in keys(ds.cubes)
        ARCEMEAnalysis.arceme_s1_position(ds, distance_threshold = 1000)
        orbit_index = ds.s1_position
        x_s1 = map(linscal, ds.time_sentinel_1_rtc.val)
        norbits = maximum(ds.s1_position)

        per_orbit = map(1:norbits) do orbit
            orbinds = findall(==(orbit), orbit_index)
            s1vals = ds.s1_indices[:, orbinds, :]
            x1, x2 = first(s1vals.time_sentinel_1_rtc) - Day(15),
            last(s1vals.time_sentinel_1_rtc) + Day(15)
            res = xmap(
                filter_spline,
                s1vals ⊘ :time_sentinel_1_rtc,
                output = XOutput(outax),
                function_args = (x_s1[orbinds], length(outax)),
            )[
                :,
                :,
                1,
                :,
            ]
            plotinds = DD.selectindices(DD.lookup(outax), x1..x2)
            res.data[1:(first(plotinds)-1), :, :] .= NaN
            res.data[(last(plotinds)+1):end, :, :] .= NaN
            res, length(plotinds) == length(outax)
        end
        fullorbits = filter(last, per_orbit)
        if isempty(fullorbits)
            # No orbit spans the full window (sparse S1): no usable S1 spline.
            nothing
        elseif length(fullorbits) > 1
            s1inds = concatenatecubes(first.(fullorbits), DD.Dim{:orbit}(1:length(fullorbits)))
            mean(s1inds, dims = :orbit)
        else
            fullorbits[1][1]
        end
    else
        # Cube has no Sentinel-1 data at all: no usable S1 spline.
        nothing
    end
    x_s2 = map(linscal, ds.time_sentinel_2_l2a.val)
    s2splines = xmap(
        filter_spline,
        ds.s2_indices ⊘ :time_sentinel_2_l2a,
        output = XOutput(outax),
        function_args = (x_s2, length(outax)),
    )[
        :,
        :,
        1,
        :,
    ]
    s1splines, s2splines
end

function expand_s2bands(cube)
    if length(cube.band_s2) == length(S2_BANDS) && all(cube.band_s2.val .== S2_BANDS)
        return cube
    else
        to_indices = [findfirst(==(v), S2_BANDS) for v in cube.band_s2.val]
        expanded = xmap(
            cube ⊘ :band_s2,
            function_args = (to_indices,),
            inplace = true,
            output = XOutput(DD.Dim{:band_s2}(collect(S2_BANDS))),
        ) do xout, xin, toinds
            xout .= missing
            xout[toinds] .= xin
        end
        permutedims(compute(expanded), ((2:ndims(expanded))..., 1))
    end
end



# S1 splines get a trailing singleton :orbit dim only when >1 full orbit exists;
# drop it when present so they are consistently [time_smooth × class × band_s1].
_squeeze_s1(s) = DD.hasdim(s, :orbit) ? dropdims(s; dims = DD.Dim{:orbit}) : s

# ---------------------------------------------------------------------------
# STEP A: per-pair, per-LC splines, veg-state label, class fraction
# ---------------------------------------------------------------------------

struct PairResult
    pair::Vector{ARCEMEAnalysis.Event}
    label::Symbol                 # :growing_season / :dormant / :skip
    s2_d::Any                     # [time_smooth × band_s2] expanded slice, or `nothing` (LC absent)
    s2_h::Any
    s1_d::Any                     # [time_smooth × band_s1] slice, or `nothing` (S1 absent)
    s1_h::Any
    fraction::Any                 # selected class fraction (metadata only)
    pair_index::Int
end

# ---------------------------------------------------------------------------
# STEP A1: per-pair spline cache (computed ONCE per pair, reused across LCs)
# ---------------------------------------------------------------------------

"""
    PairSplines

Cache of the per-pair, LC-independent spline computation.
`arceme_recovery_splines` opens both cubes and computes the dimension-agnostic
S2/S1 splines once; the land-cover loop then slices and labels from this cache
instead of recomputing. `s2_d`/`s2_h` are the S2 spline cubes expanded to the
full S2_BANDS axis (so the band axis is identical across sites, as
concatenatecubes requires). `s1_d`/`s1_h` are `nothing` when the cube has no
usable Sentinel-1 spline. `gidx` is the grouper index chosen once per pair
(symmetric across d/dhp). `ds` is kept so the LC phase can read the class axis /
fractions.
"""
struct PairSplines
    pair::Vector{ARCEMEAnalysis.Event}
    ds::Vector{Any}               # opened cubes (class axis + fractions, presence)
    s2_d::Any                     # [time_smooth × class × band_s2] expanded spline
    s2_h::Any
    s1_d::Any                     # [time_smooth × class × band_s1], or `nothing`
    s1_h::Any
    gidx::String                  # grouper index, chosen once per pair
    pair_index::Int
end

"""
    arceme_recovery_splines(pair, pair_index)

Compute the LC-independent splines for one event pair (the two `Event`s `d` and
`dhp`). Opens both cubes and returns a `PairSplines` cache. This is the per-pair
unit of `arceme_recovery_cache`; callers normally build the whole cache once and
reuse it across land-cover passes.
"""
function arceme_recovery_splines(pair, pair_index)
    ds = ARCEMEAnalysis.arceme_open.(pair)
    s1_d, s2_d = compute_indexsplines(ds[1])
    s1_h, s2_h = compute_indexsplines(ds[2])
    # Grouper index chosen symmetrically across the pair: use kNDVI only when BOTH
    # cubes have it, else NDVI (always present). This prevents the reconciliation
    # from comparing kNDVI (one cube) against NDVI (the other), which are
    # different indices with different scales. Checked on the raw (unexpanded)
    # band axis so we detect whether the source actually carried kNDVI.
    gidx = ("kNDVI" in s2_d.band_s2.val && "kNDVI" in s2_h.band_s2.val) ? "kNDVI" : "NDVI"
    # Expand each S2 spline cube to the full S2_BANDS axis so every site's slice
    # shares an identical band axis (required for concatenatecubes along :site).
    # LC-independent, so done once here rather than per LC.
    s2_d = expand_s2bands(s2_d)
    s2_h = expand_s2bands(s2_h)
    s1_d = s1_d === nothing ? nothing : _squeeze_s1(s1_d)
    s1_h = s1_h === nothing ? nothing : _squeeze_s1(s1_h)
    PairSplines(pair, ds, s2_d, s2_h, s1_d, s1_h, gidx, pair_index)
end

# ---------------------------------------------------------------------------
# STEP A2: per-pair × LC labeling (pure over the PairSplines cache)
# ---------------------------------------------------------------------------

"""
    label_pair(ps::PairSplines, selector, grouper)

Step A for one event pair and one land-cover selector, operating purely on a
cached `PairSplines` (no spline recomputation). Slices the selected LC, runs the
per-LC d/dhp veg-state reconciliation, and returns a `PairResult` whose `label`
is `:skip` when the LC is absent from a cube (discarded, `s2_d` is `nothing`) or
the reconciliation genuinely disagrees. `pair_index` identifies the site for the
`:site` axis (pairing in Step C).
"""
function label_pair(ps::PairSplines, selector, grouper)
    cls = selector.class
    if !(present(ps.ds[1], selector) && present(ps.ds[2], selector))
        # LC not present in one or both cubes: discard the site for this LC.
        return PairResult(
            ps.pair,
            :skip,
            nothing,
            nothing,
            nothing,
            nothing,
            missing,
            ps.pair_index,
        )
    end
    sl_d = select_lc(ps.s2_d, cls)
    sl_h = select_lc(ps.s2_h, cls)
    g_d = _grouper_curve(ps.s2_d, cls, ps.gidx)
    g_h = _grouper_curve(ps.s2_h, cls, ps.gidx)
    label = vegstate_group(grouper, g_d, g_h, EVENT_INDEX)
    s1d = ps.s1_d === nothing ? nothing : ps.s1_d[class=DD.At(cls)]
    s1h = ps.s1_h === nothing ? nothing : ps.s1_h[class=DD.At(cls)]
    PairResult(
        ps.pair,
        label,
        sl_d,
        sl_h,
        s1d,
        s1h,
        class_fraction(ps.ds[1], cls),
        ps.pair_index,
    )
end


# ---------------------------------------------------------------------------
# STEP A cache driver + STEP B aggregation
# ---------------------------------------------------------------------------

"""
    arceme_recovery_cache(pairs_idx)

Open every pair in `pairs_idx` and compute the LC-independent splines once,
returning a `Vector{PairSplines}`. This cache is shared across all land-cover
passes so splines are never recomputed per LC. `pairs_idx` is an iterable of
1-based indices into `arceme_validpairs()` (defaults to all valid pairs).
"""
function arceme_recovery_cache(pairs_idx = 1:length(ARCEMEAnalysis.arceme_validpairs()))
    pairs = ARCEMEAnalysis.arceme_validpairs()
    [arceme_recovery_splines(pairs[i], i) for i in pairs_idx]
end

"""
    label_all(cache, selector, grouper)

Run the per-pair × LC labeling over a cached spline store, returning a
`Vector{PairResult}`.
"""
label_all(cache, selector, grouper = GrowingSeasonSplitter()) =
    [label_pair(ps, selector, grouper) for ps in cache]

# Materialize a concatenated (DiskStack-backed) cube to an Array before xmap.
# KNOWN SEGFAULT: xmap over the DiskStack produced by concatenatecubes crashes
# reproducibly. Workaround only — see RECOVERY_PLAN.md "Known issues"; this
# needs an MWE and an upstream bug report against DiskArrayEngine/DiskArrays or
# YAXArrays.
_materialize_for_xmap(cubes_stacked) =
    YAXArray(DD.dims(cubes_stacked), Array(DD.data(cubes_stacked)))

"""
    stack_group(results; target_label)

Stack the per-site S2 spline slices of the sites whose per-LC label equals
`target_label` along a new `:pair_index` axis, relabelling the time axis to
days-since-event (`REL_TIME_AXIS`). Both the drought (`s2_d`) and drought+HP
(`s2_h`) cubes are stacked and returned **paired**: both are built from the same
`results` entries in the same order, so `D[:,:,s]` and `H[:,:,s]` are the same
site. Returns a `NamedTuple` `(D, H, sitedim)` with `D`/`H` each
`[time_smooth × band_s2 × pair_index]`. Returns `nothing` when no site matches.
`results` is the per-pair `PairResult` list for one land cover (see
`label_all`). S1 is deliberately not stacked here (deferred).
"""
function stack_group(results; target_label)
    res = filter(r -> r.label == target_label, results)
    isempty(res) && return nothing
    sitedim = DD.Dim{:pair_index}([r.pair_index for r in res])
    D = _materialize_for_xmap(
        concatenatecubes([DD.setdims(r.s2_d, REL_TIME_AXIS) for r in res], sitedim),
    )
    H = _materialize_for_xmap(
        concatenatecubes([DD.setdims(r.s2_h, REL_TIME_AXIS) for r in res], sitedim),
    )
    (D = D, H = H, sitedim = sitedim)
end

"""
    reduce_group(cube; min_samples=MIN_SAMPLES)

Reduce one stacked per-site cube (`[time_smooth × band_s2 × pair_index]`) over
the `:pair_index` axis to a cross-site summary. Returns a
`[stat × time_smooth × band_s2]` cube with `stat ∈ {"Mean","Std","SE","N"}`
(non-missing counts), leaving cells with fewer than `min_samples` valid sites
`missing`. Applied to `D` and `H` independently, so their `N` arrays may differ.
"""
function reduce_group(cube; min_samples = MIN_SAMPLES)
    # Ensure the input carries the `Missing` union so XOutput (which derives the
    # output eltype from dataeltype(ars[1])) allocates a union-typed array that
    # can hold the `missing` written for cells below `min_samples`. A caller may
    # pass a `Float64`-typed cube (e.g. a year-diff cube whose slices happened to
    # contain no `missing`), and without this the `xout .= missing` branch throws
    # a `MethodError: Cannot convert Missing to Float64`.
    data = eltype(cube.data) <: Missing ? cube.data : convert(Array{Union{Missing, Float64}}, cube.data)
    c = YAXArray(DD.dims(cube), data)
    reduced = xmap(
        c ⊘ :pair_index,
        output = XOutput(DD.Dim{:stat}(["Mean", "Std", "SE", "N"])),
    ) do xout, xin
        xclean = collect(skipmissing(xin))
        length(xclean) < min_samples && return xout .= missing
        xout .= mean(xclean),
        std(xclean),
        sem(xclean),
        length(xclean)
    end
    # xmap keeps the reduced :pair_index axis as a size-1 singleton; drop it so the
    # result is the clean [stat × time_smooth × band_s2] cross-site summary.
    Base.dropdims(reduced; dims = DD.Dim{:pair_index})
end

"""
    aggregate_group(results; target_label, transform=identity)

Top-level per-LC veg-state-group driver. For the sites whose per-LC label equals
`target_label`, stacks the paired `d` and `dhp` per-site spline slices
(`stack_group`) and reduces each over sites (`reduce_group`). Returns a
`NamedTuple` `(D, H, sitedim, d, h)`:
- `D`, `H` — unreduced paired per-site cubes `[time_smooth × band_s2 × site]`
  (the Step C seam, so the per-day paired test needs no re-stacking),
- `sitedim` — the `:pair_index` site axis (site identity),
- `d`, `h` — reduced `[stat × time_smooth × band_s2]` d and dhp curves.
Returns `nothing` when no site matches. S1 is deferred.

`transform` maps each stacked per-site cube before reduction: `transform=identity`
(the default) aggregates the raw full-window curves, while
`transform=year_diff` aggregates the year-after-minus-year-before delta (each
output cell at day `d ∈ 0:365` is `after(d) − before(d)`). A single function
replaces the former raw-only and year-diff-only aggregation helpers.
"""
function aggregate_group(results; target_label, transform = identity)
    st = stack_group(results; target_label)
    st === nothing && return nothing
    D = transform(st.D)
    H = transform(st.H)
    (; D = D, H = H, sitedim = st.sitedim, d = reduce_group(D), h = reduce_group(H))
end

"""
    year_diff(cube)

Map a stacked per-site cube (`[time_smooth × band_s2 × pair_index]`, full
`-365:365` axis) to the year-after-minus-year-before difference: each output cell
at day `d ∈ 0:365` is `after(d) − before(d)`, pairing day `d` (after) with day
`d − 365` (before), i.e. always 365 days apart. `missing` is preserved pointwise
(a site/band where either half is missing yields `missing`). Returns a cube with
the `time_smooth` axis relabelled to `0:365` (366 points). Applied to the stacked
`D`/`H` cubes, so the per-site pairing needed by Step C survives.
"""
function year_diff(cube)
    ts = cube.time_smooth
    dayvals = collect(ts.val)
    # After half = days 0:365, before half = days -365:0 (day 0 appears in both,
    # so both halves have 366 points); each after point pairs 365 days on.
    iafter = findall(≥(0), dayvals)
    ibefore = findall(≤(0), dayvals)
    @assert length(iafter) == length(ibefore)
    newts = DD.Dim{:time_smooth}(0:(length(iafter)-1))
    # Slice-subtract along time, keep band × site axes.
    after = cube[time_smooth=DD.At(dayvals[iafter])]
    before = cube[time_smooth=DD.At(dayvals[ibefore])]
    diffdata = Array(after.data) .- Array(before.data)
    YAXArray((newts, cube.band_s2, cube.pair_index), diffdata)
end

# ---------------------------------------------------------------------------
# STEP C: per-day paired significance of d - dhp on the recovery window
# ---------------------------------------------------------------------------

# One-sided p is not used; the two-sided signed-rank tests whether the d/dhp
# difference at a (day, band) point is centered on zero across sites. Per-point
# missing policy (decision 3): a site contributes only where BOTH d and dhp are
# non-missing at that exact day. Fewer than 2 usable sites -> `missing` (the
# empty/`<2` case is untestable; HypothesisTests.SignedRankTest errors on empty
# input). HypothesisTests handles all-equal differences correctly (e.g. n=5 all
# one sign -> exact p = 2/2^5), so no separate all-ties guard is needed.
function _signed_rank_p(diffs)
    usable = collect(skipmissing(diffs))
    length(usable) < 2 && return missing
    HypothesisTests.pvalue(HypothesisTests.SignedRankTest(usable))
end

"""
    recovery_pvalues(st; bands=S2_BANDS)

Raw per-day, per-band paired significance of d - dhp on the recovery window. The
recovery window is derived from `st.D`'s `time_smooth` axis as the positions where
the day value is `≥ 0` (so it works for both the full `-365:365` axis and the
`year_diff` `0:365` axis). For each (day, band), forms the site-level difference
`D[day,band,:] .- H[day,band,:]` (pointwise-missing per decision 3) and runs a
two-sided Wilcoxon signed-rank test across sites. `st` is an `aggregate_group`
output. Returns `(p, n, days, bands)`: `p` is a `[ntest × band_s2]` array of raw
p-values (`missing` where untestable), `n` the `[ntest × band_s2]` count of usable
sites per point, and `days` the recovery-window day values (e.g. `0:365`).
"""
function recovery_pvalues(st; bands = S2_BANDS)
    dayvals = collect(st.D.time_smooth.val)
    days = findall(≥(0), dayvals)          # positions on the axis (recovery window)
    ibands = [findfirst(==(b), st.D.band_s2.val) for b in bands]
    D = Array(st.D.data)[:, ibands, :]
    H = Array(st.H.data)[:, ibands, :]
    Dsub = D[days, :, :]
    Hsub = H[days, :, :]
    ndays = length(days)
    nband = length(bands)
    p = Array{Union{Missing,Float64}}(missing, ndays, nband)
    n = zeros(Int, ndays, nband)
    for j = 1:nband, i = 1:ndays
        diff = Dsub[i, j, :] .- Hsub[i, j, :]
        usable = collect(skipmissing(diff))
        n[i, j] = length(usable)
        p[i, j] = _signed_rank_p(diff)
    end
    (p = p, n = n, days = dayvals[days], bands = bands)
end

"""
    fdr_bh(pvals; α=FDR_ALPHA)

Benjamini-Hochberg false-discovery-rate control over the flat list of raw
p-values (`missing`/untestable points excluded). Returns `(q, significant)` with
`q` the adjusted values (same shape as `pvals`, `missing` preserved) and
`significant` a same-shape `Bool`/`missing` mask, true where `q ≤ α`.
"""
function fdr_bh(pvals; α = FDR_ALPHA)
    flat = vec(pvals)
    finite_inds = findall(!ismissing, flat)
    n = length(finite_inds)
    q = similar(flat, Union{Missing,Float64})
    sig = similar(flat, Union{Missing,Bool})
    q .= missing
    sig .= missing
    if n > 0
        ps = Float64.(flat[finite_inds])
        order = sortperm(ps)                # rank (smallest first) -> position in ps
        sorted = ps[order]
        m = sorted .* n ./ (1:n)            # raw BH multipliers p_(i) * n / i
        qs = similar(m)
        qs[n] = m[n]
        for r = (n-1):-1:1               # backward cumulative min (step-up)
            qs[r] = min(m[r], qs[r+1])
        end
        for (rank, k) in enumerate(order)
            flatidx = finite_inds[k]
            q[flatidx] = qs[rank]
            sig[flatidx] = qs[rank] ≤ α
        end
    end
    (q = reshape(q, size(pvals)), significant = reshape(sig, size(pvals)))
end

"""
    arceme_recovery_test(st; bands=S2_BANDS, α=FDR_ALPHA)

Top-level Step C driver over an `aggregate_group` output `st`. Computes
the raw per-day per-band signed-rank p-values on the recovery window, applies
Benjamini-Hochberg FDR across all (day, band) points, and returns a `NamedTuple`
`(p, q, significant, n, days, bands)` for Step D's significance shading.
"""
function arceme_recovery_test(st; bands = S2_BANDS, α = FDR_ALPHA)
    raw = recovery_pvalues(st; bands)
    bh = fdr_bh(raw.p; α)
    (;
        p = raw.p,
        q = bh.q,
        significant = bh.significant,
        n = raw.n,
        days = raw.days,
        bands = raw.bands,
    )
end

"""
    arceme_recovery_pairs(cache, lc, grouper, target_label)

Per-(land-cover, veg-state-group) aggregation, the single place that computes the
raw + year-diff aggregates (shared by the full-run driver so the aggregation logic
is not duplicated). Labels every pair in `cache` for `lc` via `grouper`, then
aggregates the sites whose label equals `target_label`: the year-after-minus-year-before
difference (`aggregate_group(...; transform=year_diff)`, exposed as `curves`/`test`)
and the raw full-window curves (`aggregate_group(...; transform=identity)`, exposed
as `curves_full`/`test_full`), each with its Step C significance
(`arceme_recovery_test`). Returns `nothing` when no site matches the group, else a
NamedTuple with `curves`/`test` (delta), `curves_full`/`test_full` (raw) and the
`target_label` site count `sites`.
"""
function arceme_recovery_pairs(cache, lc, grouper, target_label)
    res = label_all(cache, LandCoverSelector(class = lc), grouper)
    grp = aggregate_group(res; target_label)
    grp === nothing && return nothing
    grpd = aggregate_group(res; target_label, transform = year_diff)
    (;
        curves = grpd,
        test = arceme_recovery_test(grpd),
        curves_full = grp,
        test_full = arceme_recovery_test(grp),
        sites = count(r -> r.label == target_label, res),
    )
end

# ---------------------------------------------------------------------------
# STEP D: plots + CSV of the reduced recovery curves and Step C significance
# ---------------------------------------------------------------------------
#
# Step C produces, per (LC, veg-state group), the reduced d/dhp mean+SE curves
# (`curves.d` / `curves.h`, [stat × time_smooth × band]) and per-day per-band
# significance (`test.p`/`test.q`/`test.significant`, [366 × 8] on the recovery
# window). Step D renders these as Makie delta panels and writes the
# underlying numbers as long/tidy CSV. Plotting uses `Makie` directly; the
# caller must load a backend (e.g. `CairoMakie` or `GLMakie`) first.

# Sanitize a land-cover name into a filename-safe token (spaces and slashes are
# not safe in file names).
_safe_token(s) = replace(strip(s), r"[/ ]+" => "_")

# The raw full-window analysis from an aggregation entry: `curves_full`/`test_full`
# when the entry is a combined raw+delta entry (as `arceme_recovery_pairs`
# produces), falling back to `curves`/`test` for a plain raw-only entry.
_raw_curves(ent) = hasproperty(ent, :curves_full) ? ent.curves_full : ent.curves
_raw_test(ent) = hasproperty(ent, :test_full) ? ent.test_full : ent.test

# Prefix each tidy CSV row with its land-cover label so the per-LC files stay
# self-describing (exactly one LC per file, but the column is retained).
_with_lc(lc, rows) = (merge((lc = lc,), r) for r in rows)

# Pull the per-band x (days-since-event), mean and SE 1-D vectors out of one
# reduced [stat × time_smooth × band] cube.
function _band_curve(red, band)
    ib = findfirst(==(band), red.band_s2.val)
    (;
        x = collect(red.time_smooth.val),
        mean_ = vec(red.data[1, :, ib]),
        se = vec(red.data[3, :, ib]),
    )
end

"""
    recovery_panel!(ax, grp, test, band; xwin=:full, yzero=false)

Draw one band's panel onto `ax` from a Step B/C output pair: the `d` (red) and
`dhp` (blue) mean curves with `Mean ± SE` ribbons, the day-0 (event) vline, and
Step C significance shading on the recovery window. `grp` is an
`aggregate_group` output, `test` an `arceme_recovery_test` output.

Keyword options:
- `xwin=:full` plots the whole `-365..365` window with the recovery window
  (0..365) shaded as context; `xwin=:recovery` zooms the x-axis to `0..365`.
- `yzero=true` centers the y-axis on zero (used for the year-diff delta panels).
"""
function recovery_panel!(ax, grp, test, band; xwin = :full, yzero = false)
    d = _band_curve(grp.d, band)
    h = _band_curve(grp.h, band)
    dayvals = test.days
    ib = findfirst(==(band), test.bands)
    p = vec(test.p[:, ib])
    sig = vec(test.significant[:, ib])

    # Plotting a vector of `missing` (a band or side with no data) crashes Makie's
    # float_type, so each side is drawn only when it has finite values. A band with
    # no data on either side (e.g. a veg-state group with fewer than MIN_SAMPLES
    # contributing sites) is left as an empty panel.
    hasval(v) = any(x -> !ismissing(x) && isfinite(x), v)
    d_has = hasval(d.mean_) && hasval(d.se)
    h_has = hasval(h.mean_) && hasval(h.se)
    (d_has || h_has) || return ax

    # Recovery window context (only shown when not already the full x-range).
    if xwin == :full
        vspan!(ax, 0, 365; color = (:grey, 0.06))
    end

    # Significance shading, per day on the recovery window.
    for (k, day) in enumerate(dayvals)
        !ismissing(p[k]) &&
            p[k] < 0.05 &&
            vspan!(ax, day, day + 1; color = (:orange, 0.15))
        !ismissing(sig[k]) && sig[k] && vspan!(ax, day, day + 1; color = (:red, 0.35))
    end

    # Mean curves + Mean ± SE ribbons, per side that has data.
    if d_has
        band!(ax, d.x, d.mean_ .- d.se, d.mean_ .+ d.se; color = (:red, 0.18))
        lines!(ax, d.x, d.mean_; color = :red, label = "d")
    end
    if h_has
        band!(ax, h.x, h.mean_ .- h.se, h.mean_ .+ h.se; color = (:blue, 0.18))
        lines!(ax, h.x, h.mean_; color = :blue, label = "dhp")
    end

    vlines!(ax, [0]; color = :black, linestyle = :dash)
    ax.title = string(band)
    xwin == :recovery && xlims!(ax, 0, 365)
    # Center the y-axis on zero so the delta curve's positive/negative excursions
    # read symmetrically about the no-difference line. The range is taken from the
    # plotted mean ± SE data, so the axis is guaranteed to include zero.
    if yzero
        finite = Float64[]
        d_has && append!(finite, filter(x -> !ismissing(x) && isfinite(x), d.mean_ .- d.se))
        d_has && append!(finite, filter(x -> !ismissing(x) && isfinite(x), d.mean_ .+ d.se))
        h_has && append!(finite, filter(x -> !ismissing(x) && isfinite(x), h.mean_ .- h.se))
        h_has && append!(finite, filter(x -> !ismissing(x) && isfinite(x), h.mean_ .+ h.se))
        m = isempty(finite) ? 1.0 : maximum(abs, finite)
        ylims!(ax, -m, m)
    end
    ax
end

"""
    _csv_rows(grp, test)

Build the long/tidy row vectors for one (LC-agnostic) analysis entry: curve rows
(`lc × scenario(d|h) × stat × day × band × value`) and test rows
(`lc × day × band × n × p × q × significant`). Shared by
`arceme_recovery_export_csv` so the raw and `year_diff` analyses reuse the same
row construction.
"""
function _csv_rows(grp, test)
    curves_rows = NamedTuple[]
    tests_rows = NamedTuple[]
    for scenario in (:d, :h)
        red = grp[scenario]
        for (bi, band) in enumerate(red.band_s2.val)
            for (ti, day) in enumerate(red.time_smooth.val)
                for (si, stat) in enumerate(["Mean", "Std", "SE", "N"])
                    push!(
                        curves_rows,
                        (
                            scenario = scenario,
                            stat = stat,
                            day = day,
                            band = band,
                            value = red.data[si, ti, bi],
                        ),
                    )
                end
            end
        end
    end
    for (i, day) in enumerate(test.days)
        for (bi, band) in enumerate(test.bands)
            push!(
                tests_rows,
                (
                    day = day,
                    band = band,
                    n = test.n[i, bi],
                    p = test.p[i, bi],
                    q = test.q[i, bi],
                    significant = test.significant[i, bi],
                ),
            )
        end
    end
    (curves_rows, tests_rows)
end

"""
    arceme_recovery_export_csv(ent; lc, splitter, group, outdir)

Write one set of Step B/C numbers as long/tidy CSV for a single land cover,
reusing `_csv_rows` for the rows. Files go under
`outdir/<splitter>/<group>/` with the `-s2` sensor token and the land-cover token
in the name, e.g. `Tree_cover-s2-recovery_curves.csv`. The raw full-window
curves/tests (from `curves_full`/`test_full`, or `curves`/`test` for a raw-only
entry) are written as `recovery_curves.csv`/`recovery_tests.csv`; a combined
raw+delta entry additionally writes the year-diff numbers (from `curves`/`test`)
as `recovery_year_diff_curves.csv`/`recovery_year_diff_tests.csv` and the full
raw numbers as `recovery_full_curves.csv`/`recovery_full_tests.csv`. Each file
retains the `lc` column (one LC per file). Returns the list of written file paths.
"""
function arceme_recovery_export_csv(ent; lc, splitter, group, outdir)
    out = joinpath(outdir, string(splitter), group_token(splitter, group))
    mkpath(out)
    base = joinpath(out, "$(_safe_token(lc))-$(_sensor_token())")

    c, t = _csv_rows(_raw_curves(ent), _raw_test(ent))
    files = String[]
    # CSV.write accepts any Tables.jl source; materialize the single-pass
    # _with_lc generator so the rows can be streamed safely.
    push!(files, CSV.write(joinpath(base * "-recovery_curves.csv"), collect(_with_lc(lc, c))))
    push!(files, CSV.write(joinpath(base * "-recovery_tests.csv"), collect(_with_lc(lc, t))))
    if hasproperty(ent, :curves_full)
        yd_c, yd_t = _csv_rows(ent.curves, ent.test)
        push!(files, CSV.write(joinpath(base * "-recovery_year_diff_curves.csv"), collect(_with_lc(lc, yd_c))))
        push!(files, CSV.write(joinpath(base * "-recovery_year_diff_tests.csv"), collect(_with_lc(lc, yd_t))))
        fl_c, fl_t = _csv_rows(ent.curves_full, ent.test_full)
        push!(files, CSV.write(joinpath(base * "-recovery_full_curves.csv"), collect(_with_lc(lc, fl_c))))
        push!(files, CSV.write(joinpath(base * "-recovery_full_tests.csv"), collect(_with_lc(lc, fl_t))))
    end
    files
end

"""
    arceme_recovery_plot_delta(ent; lc, splitter, group, outdir)

Combined figure for a `diff_year=true` entry: one row per band, with the
full-window raw curves (`curves_full`, `xwin=:full`) in the left column and the
year-after-minus-year-before delta curve (`curves`, `xwin=:recovery`) in the
right column, so the seasonal-subtracted delta can be read against its full
context. Saves `{lc}-s2-recovery_delta.png` under `outdir/<splitter>/<group>/`
and returns its path. `ent` is a combined raw+delta entry, `lc` the land-cover
label, `splitter` the `VegStateGrouper`, `group` the veg-state group `Symbol`,
and `outdir` the root plots directory.
"""
function arceme_recovery_plot_delta(ent; lc, splitter, group, outdir)
    grp = ent.curves
    test = ent.test
    grpf = ent.curves_full
    testf = ent.test_full
    bands = collect(grp.d.band_s2.val)
    nb = length(bands)
    fig = Figure(size = (640, 170 * nb))
    for (bi, band) in enumerate(bands)
        ax_full = Axis(fig[bi, 1], title = string(band))
        recovery_panel!(ax_full, grpf, testf, band; xwin = :full)
        ax_delta = Axis(fig[bi, 2])
        recovery_panel!(ax_delta, grp, test, band; xwin = :recovery, yzero = true)
    end
    Label(fig[0, 1], "Full window", tellwidth = false)
    Label(fig[0, 2], "Year-after minus year-before", tellwidth = false)
    out = joinpath(outdir, string(splitter), group_token(splitter, group))
    mkpath(out)
    fname = joinpath(out, "$(_safe_token(lc))-$(_sensor_token())-recovery_delta.png")
    save(fname, fig)
    fname
end

"""
    arceme_recovery_plot(ent; lc, splitter, group, outdir, data_outdir=outdir)

Top-level Step D driver for one (LC, veg-state group) entry: render the
full-vs-delta figure (`arceme_recovery_plot_delta`) for a combined raw+delta
entry, and write the per-LC CSV tables (`arceme_recovery_export_csv`). Figures
go under `outdir/<splitter>/<group>/`, CSV under
`data_outdir/<splitter>/<group>/`. Returns `(figures, data)` — the lists of
written figure and CSV file paths.
"""
function arceme_recovery_plot(ent; lc, splitter, group, outdir, data_outdir = outdir)
    figures = String[]
    if hasproperty(ent, :curves_full)
        push!(figures, arceme_recovery_plot_delta(ent; lc, splitter, group, outdir))
    end
    data = arceme_recovery_export_csv(ent; lc, splitter, group, outdir = data_outdir)
    (figures = figures, data = data)
end

const FULL_LC_SET = ["Tree cover", "Shrubland", "Grassland", "Cropland", "Bare /sparse vegetation"]

"""
    arceme_run_recovery(cache; lc_set=FULL_LC_SET, grouper=GrowingSeasonSplitter(),
                        groups=vegstate_groups(grouper), plots_root=plots_root,
                        data_root=data_root)

Unified artifact-producing runner. It builds no cache itself: the caller passes a
pre-built spline `cache` (see `arceme_recovery_cache`) covering exactly the range
of pairs to analyze. Loops over each land-cover class in `lc_set` and each
veg-state group in `groups` (defaults to the grouper's declared groups, see
`vegstate_groups`). For each (LC, group) the raw + year-diff aggregation is
delegated to `arceme_recovery_pairs` (the single place that computes both
aggregates); groups with no sites for an LC are skipped with a report note. For
each qualifying (LC, group) `arceme_recovery_plot` writes the delta figure into
`plots_root/<splitter>/<group>/` and the per-LC CSV tables into
`data_root/<splitter>/<group>/`.

Returns a `report` Dict mapping each (LC, group) key to `(sites, figures, data)`
— the veg-state-group site count and the lists of written figure / data file
paths — so the run can be audited.
"""
function arceme_run_recovery(cache;
    lc_set = FULL_LC_SET,
    grouper = GrowingSeasonSplitter(),
    groups = vegstate_groups(grouper),
    plots_root = joinpath(pwd(), "recovery_plots"),
    data_root = joinpath(pwd(), "recovery_data"),
)
    report = Dict{String,Any}()
    for lc in lc_set
        for group in groups
            key = "$(_safe_token(lc)) / $(group_token(grouper, group))"
            agg = arceme_recovery_pairs(cache, lc, grouper, group)
            if agg === nothing
                report[key] = (sites = 0, figures = String[], data = String[])
                continue
            end
            figdata = arceme_recovery_plot(
                agg;
                lc,
                splitter = grouper,
                group,
                outdir = plots_root,
                data_outdir = data_root,
            )
            report[key] = (sites = agg.sites, figures = figdata.figures, data = figdata.data)
        end
    end
    report
end

end # module
