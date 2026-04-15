export arceme_bias_corrected_fp, arceme_uncorrected_fp

function sincosparams(c,timname = :time_sentinel_2_l2a)
    fitmat = bare_matrix(c,timname)
    xmap(
        fit_sincos, 
        c ⊘ timname, 
        output=XOutput(DD.Dim{:param}(["p1","p2","p3"]),outtype=Float32),
        inplace=false,
        function_args=(fitmat,),
        allow_threads=true
    )
end

function fit_sincos(ts,fitmat)
    igood = findall(!isnan,ts)
    view(fitmat,igood,:) \ view(ts,igood)
end

function bare_matrix(c,timname)
    ms_per_year = 365.25*24*60*60*1000
    timdim = DD.dims(c,timname)
    t = map(timdim.val) do x
        (x-first(timdim)).value/ms_per_year*2pi
    end
    [sin.(t) cos.(t) fill(1.0,length(t))]
end

function nstrata(strata)
    strata=="ESA_LC" && return length(arceme_classes)
    strata == "CTY" && return 18
    strata == "MCTY" && return length(HRL.hrl_legends["MCTY"])
end

function classaxis(strata, ncl)
    if strata == "ESA_LC"
        DD.Dim{:class}(collect(values(arceme_classes))[1:ncl])
    else
        DD.Dim{:class}(collect(values(HRL.hrl_legends[strata]))[1:ncl])
    end
end


"""
`arceme_bias_corrected_fp(band, dataset; lccube=lckeymap.(dataset.ESA_LC[time=1]), ncl=12, timeaxis=:time_sentinel_2_l2a)`

Computes a cloud-biased corrected footprint of `dataset[band]` aggregated by stratification class for the 
provided stratification cube `lccube` (`ncl` classes). 
"""
function arceme_bias_corrected_fp(band::String, dataset::Dataset; strata="ESA_LC", timeaxis=:time_sentinel_2_l2a)
    ncl=nstrata(strata)
    lccube = lckeymap(dataset, strata=strata)
    cloudcube = dataset.cloud_mask
    sclcube = dataset.SCL
    inputcube = dataset[band]
    pars = ARCEMEAnalysis.sincosparams(inputcube)[:,:,:,1]
    timdim = DD.dims(inputcube,timeaxis)
    fitmat = YAXArray((timdim, pars.param),ARCEMEAnalysis.bare_matrix(inputcube,timeaxis))

    # clearsky_expected = xmap(pars ⊘ :param, fitmat ⊘ :param, inplace=false, lazy=true) do p, t
    #     p[1]*t[1]+p[2]*t[2]+p[3]*t[3]
    # end
    # win_clearsky = YAXArrays.windows(clearsky_expected, lccube, expected_groups=1:ncl)
    # fp_clearsky_expected = YAXArrays.compute(mean.(win_clearsky).data)[:, 1, 1, 1, :]

    # clouded_expected = xmap(pars ⊘ :param, fitmat ⊘ :param, cloudcube, sclcube, inplace=false) do p, t, cl, scl
    #     if _is_cloud(cl,scl)
    #         return NaN
    #     else
    #         p[1]*t[1]+p[2]*t[2]+p[3]*t[3]
    #     end
    # end
    # win_clouded = YAXArrays.windows(clouded_expected, lccube, expected_groups=1:ncl)
    # fp_clouded_expected = YAXArrays.compute(mean.(win_clouded).data)[:, 1, 1, 1, :]
    varax = DD.Dim{:Variable}(["expected", "expected_cloudfiltered", "data_cloudfiltered"])
    predcube = xmap(pars ⊘ :param, fitmat ⊘ :param, inputcube, cloudcube, sclcube, inplace=true, output=XOutput(varax, outtype=Float32)) do out, p, t, x, cl, scl
        expected = p[1] * t[1] + p[2] * t[2] + p[3] * t[3]
        out[1] = expected
        if _is_cloud(cl,scl)
            out[2] = NaN
            out[3] = NaN
        else
            out[2] = expected
            out[3] = x
        end
    end

    win_data = YAXArrays.windows(predcube, lccube, expected_groups=1:ncl)
    fp = YAXArrays.compute(mean.(win_data).data, showprogress=false)[:, :, 1, 1, 1, :]
    # win_data = YAXArrays.windows(inputcube_filtered, lccube, expected_groups=1:ncl)
    # fp = YAXArrays.compute(mean.(win_data).data)[:, 1, 1, :]

    newdata = fp[:, 3, :] .+ fp[:, 1, :] .- fp[:, 2, :]

    classax = classaxis(strata, ncl)

    Dataset(
        fp = YAXArray((classax,timdim),newdata),
        fp_uncorrected=YAXArray((classax, timdim), fp[:, 3, :]),
        fp_clearsky_expected=YAXArray((classax, timdim), fp[:, 1, :]),
        fp_clouded_expected=YAXArray((classax, timdim), fp[:, 2, :]),
        params=pars,
        smooth_matrix=fitmat,
    )
    
end

"""
`arceme_uncorrected_fp(band, dataset; strata="ESA_LC", ncl=nstrata(strata), timeaxis=:time_sentinel_1_rtc)`

Computes a cloud-biased corrected footprint aggregated by stratification class for the 
provided inputcube, cloud mask and stratification cube (`ncl` classes). 
"""
function arceme_uncorrected_fp(band, dataset; strata="ESA_LC", timeaxis=:time_sentinel_1_rtc)
    ncl=nstrata(strata)
    lccube = lckeymap(dataset, strata=strata)
    inputcube = dataset[band]
    timdim = DD.dims(inputcube,timeaxis)
    
    win = YAXArrays.windows(inputcube, lccube, expected_groups=1:ncl)
    fp = mean.(win)[:,1,1,:].data

    abundance = counter(lccube)
    classax = classaxis(strata, ncl)


    Dataset(
        fp=YAXArray((classax, timdim), fp),
        class_fractions=YAXArray((classax,), [abundance[i] for i in 1:ncl])
    )
    
end