import SpectralIndices as SI
normalize_s2name(n) = match(r"B\d",n) === nothing ? Symbol(n) : Symbol("B0$(last(n))")
struct NTWrapper{F,names,C,I} <: Function
    f::F
    names::Val{names}
    consts::C
    indices::I
end
function (ntw::NTWrapper{F,names})(cl,scl,args...) where {F,names} 
    ntw.f(cl,scl,NamedTuple{names}(args),ntw.consts,ntw.indices) 
end
NTWrapper(f,names::Tuple,consts::NamedTuple,indices) = NTWrapper(f,Val(names),consts,indices)
function inner_compute_indices(cl,scl,bands,consts,indices_tuple)
    (cl > 0 || (scl in (1, 3, 7, 8, 9, 10, 11))) && return map(_ -> NaN, indices_tuple)
    bandparams = map(ARCEMEAnalysis.boa,bands)
    allparams = (; bandparams..., consts...)
    listofindices(indices_tuple, allparams)
end

function _compute_indices(ds,indices,platform)
    indices_tuple = ((SI.indices[k] for k in indices)...,)

    allvars = mapreduce(r->Set(string.(SI._band_names(r))),union!,indices_tuple)

    needed_bands = intersect(allvars,keys(SpectralIndices.bands))
    needed_constants = intersect(allvars,keys(SpectralIndices.constants))
    undefined_constants = setdiff(allvars,needed_bands, needed_constants)
    default_values = (;(Symbol(c)=>SI.constants[c].default for c in needed_constants)...)
    bands_to_load = (sort([normalize_s2name(SI.bands[b].platforms[platform].band) for b in needed_bands])...,)

    wrapped_function = NTWrapper(inner_compute_indices, (Symbol.(needed_bands)...,),default_values,indices_tuple)

    bandargs = map(b->ds.cubes[b],bands_to_load)
    output = map(_ -> XOutput(; outtype=Float32), indices_tuple)


    return NamedTuple{(Symbol.(indices)...,)}(xmap(wrapped_function, ds.cloud_mask, ds.SCL, bandargs...; output, inplace=false))
end