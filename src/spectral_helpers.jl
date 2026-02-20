import SpectralIndices as SI
import DiskArrayEngine as DAE
import OnlineStats
import CSV

function load_sigmadata()
    CSV.read(joinpath(@__DIR__, "..", "data", "sigmalist.csv"), NamedTuple)
end

const sigmadata = load_sigmadata()


normalize_name(n) = match(r"B\d", n) !== nothing ? Symbol("B0$(last(n))") :
                    in(n, ("VH", "VV")) ? Symbol(lowercase(n)) : Symbol(n)
struct NTWrapper{F,names,C,I} <: Function
    f::F
    names::Val{names}
    consts::C
    indices::I
end
function (ntw::NTWrapper{F,names})(args...) where {F,names}
    ntw.f(NamedTuple{names}(args), ntw.consts, ntw.indices)
end
NTWrapper(f,names::Tuple,consts::NamedTuple,indices) = NTWrapper(f,Val(names),consts,indices)
function inner_compute_indices_s2(bands, consts, indices_tuple)
    (bands.cl > 0 || (bands.scl in (1, 3, 7, 8, 9, 10, 11))) && return map(_ -> NaN, indices_tuple)
    bandparams = map(ARCEMEAnalysis.boa,bands)
    allparams = (; bandparams..., consts...)
    listofindices(indices_tuple, allparams)
end
function inner_compute_indices_s1(bands, consts, indices_tuple)
    allparams = (; bands..., consts...)
    listofindices(indices_tuple, allparams)
end

function _compute_indices(ds,indices,platform)
    indices_tuple = ((SI.indices[k] for k in indices)...,)

    allvars = mapreduce(r->Set(string.(SI._band_names(r))),union!,indices_tuple)

    needed_bands = intersect(allvars,keys(SpectralIndices.bands))
    needed_constants = intersect(allvars,keys(SpectralIndices.constants))
    undefined_constants = setdiff(allvars,needed_bands, needed_constants)
    default_values = (;(Symbol(c)=>SI.constants[c].default for c in needed_constants)...)
    bands_to_load = (sort([normalize_name(SI.bands[b].platforms[platform].band) for b in needed_bands])...,)
    bandargs = map(b -> ds.cubes[b], bands_to_load)
    output = map(_ -> XOutput(; outtype=Float32), indices_tuple)
    if startswith(platform, "sentinel2")
        wrapped_function = NTWrapper(inner_compute_indices_s2, (:cl, :scl, Symbol.(needed_bands)...), default_values, indices_tuple)
        NamedTuple{(Symbol.(indices)...,)}(xmap(wrapped_function, ds.cloud_mask, ds.SCL, bandargs...; output, inplace=false))
    elseif startswith(platform, "sentinel1")
        wrapped_function = NTWrapper(inner_compute_indices_s1, (Symbol.(needed_bands)...,), default_values, indices_tuple)
        res = xmap(wrapped_function, bandargs...; output, inplace=false)
        length(indices) == 1 && (res = (res,))
        NamedTuple{(Symbol.(indices)...,)}(res)
    end
end

#Helper functions to compute the indices from named tuples for type stability
compute_indexx(index::SpectralIndices.SpectralIndex{<:Any,B}, params::NamedTuple) where B = index.compute(Float64, params[B]...)
function listofindices(indices, values)
    map(indices) do index
        compute_indexx(index, values)
    end
end

#Inner function to compute absolute difference of N and R band
function absnormdiff(b1, b2, cl, scl)
    if _is_cloud(cl, scl)
        NaN
    else
        abs(boa(b1) - boa(b2))
    end
end


function estimate_sigma_per_lc(ds;bands=(:B08,:B04))
    stat = DAE.DerivedOnlineStat{OnlineStats.ExpandingHist,OnlineStats.median,OnlineStats.fit!,(200,)}
    f = DAE.disk_onlinestat(stat, Union{Float64,Missing}, identity, ARCEMEAnalysis.lckeymap)

    si = size(ds[first(bands)].data)
    iw1 = DAE.InputArray(absnormdiff.(ds[first(bands)].data, ds[last(bands)].data, ds.cloud_mask.data, ds.SCL.data))
    iw2 = DAE.InputArray(ds.ESA_LC[time=1].data)
    ow = DAE.create_outwindows((si..., 12), windows=(fill.(1, si)..., [1:12]))
    op = DAE.GMDWop((iw1, iw2), (ow,), f)
    res = DAE.compute(DAE.results_as_diskarrays(op)[1])
    res[1, 1, 1, :]
end


function _kNDVI(cl, scl, nir, red, lc, sigmadata)
    _is_cloud(cl, scl) && return NaN
    sigma = sigmadata[lckeymap(lc)]
    banddiff = boa(nir) - boa(red)
    return tanh(((banddiff) / (2 * sigma))^2) * sign(banddiff)
end

"""

"""
function arceme_kndvi(ds, sigmadata=estimate_sigma_per_lc(ds))
    kndvi = xmap(_kNDVI, ds.cloud_mask, ds.SCL, ds.B08, ds.B04, ds.ESA_LC[time=1], output=XOutput(outtype=Float32), function_args=(sigmadata,), inplace=false)
    ds.cubes[:kndvi] = kndvi
end



function arceme_kndvi_pair(ds1, ds2)
    missmean(x1, x2) =
        if ismissing(x1)
            ismissing(x2) ? missing : x2
        else
            ismissing(x2) ? x1 : 0.5 * (x1 + x2)
        end
    sigmadata1 = estimate_sigma_per_lc(ds1)
    sigmadata2 = estimate_sigma_per_lc(ds2)
    sigmadata = missmean.(sigmadata1, sigmadata2)
    arceme_kndvi(ds1, sigmadata)
    arceme_kndvi(ds2, sigmadata)
end
