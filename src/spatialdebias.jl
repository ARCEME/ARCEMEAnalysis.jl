export arceme_bias_corrected_fp

function sincosparams(c,timname = :time_sentinel_2_l2a)
    fitmat = bare_matrix(c,timname)
    xmap(
        fit_sincos, 
        c ⊘ timname, 
        output=XOutput(DD.Dim{:param}(["p1","p2","p3"]),outtype=Float32),
        inplace=false,
        function_args=(fitmat,)
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

"""
    arceme_bias_corrected_fp(inputcube,cloudcube,lccube)

Computes a cloud-biased corrected footprint aggregated by land cover class for the 
provided inputcube, cloud mask and land cover cube. 
"""
function arceme_bias_corrected_fp(inputcube, dataset)
    lccube = dataset.ESA_LC[time=1]
    cloudcube = dataset.cloud_mask
    sclcube = dataset.SCL
    pars = ARCEMEAnalysis.sincosparams(inputcube)[:,:,:,1]

    fitmat = YAXArray((inputcube.time_sentinel_2_l2a, pars.param),ARCEMEAnalysis.bare_matrix(inputcube,:time_sentinel_2_l2a))

    clearsky_expected = xmap(pars ⊘ :param, fitmat ⊘ :param,inplace=false) do p,t
        p[1]*t[1]+p[2]*t[2]+p[3]*t[3]
    end
    win_clearsky = YAXArrays.windows(clearsky_expected,lccube,expected_groups=0:100)
    fp_clearsky_expected = mean.(win_clearsky)[:,1,1,1,:].data

    clouded_expected = xmap(pars ⊘ :param, fitmat ⊘ :param, cloudcube, sclcube, inplace=false) do p, t, cl, scl
        if (cl > 0 || scl in (1, 3, 8, 9, 10, 11))
            return NaN
        else
            p[1]*t[1]+p[2]*t[2]+p[3]*t[3]
        end
    end
    win_clouded = YAXArrays.windows(clouded_expected,lccube,expected_groups=0:100)
    fp_clouded_expected = mean.(win_clouded)[:,1,1,1,:].data

    inputcube_filtered = xmap(inputcube, cloudcube,sclcube,inplace=false,output=XOutput(outtype=Float32)) do x,cl,scl
        (cl > 0 || scl in (1, 3, 8, 9, 10, 11)) ? NaN : x
    end
    win_data = YAXArrays.windows(inputcube_filtered,lccube,expected_groups=0:100)
    fp = mean.(win_data)[:,1,1,:].data

    classkeys = collect(keys(arceme_classes))[2:end]
    abundance = counter(lccube)
    sort!(classkeys,by=i->(abundance[i],i),rev=true)
    newdata = fp[classkeys,:] .+ fp_clearsky_expected[classkeys,:] .- fp_clouded_expected[classkeys,:]
    classax = DD.Dim{:lc}([arceme_classes[i] for i in classkeys])

    Dataset(
        fp = YAXArray((classax,inputcube.time_sentinel_2_l2a),newdata),
        fp_uncorrected=YAXArray((classax, inputcube.time_sentinel_2_l2a), fp[classkeys, :]),
        clearsky_expected = YAXArray((classax,inputcube.time_sentinel_2_l2a),fp_clearsky_expected[classkeys,:]),
        clouded_expected = YAXArray((classax,inputcube.time_sentinel_2_l2a),fp_clouded_expected[classkeys,:]),
        lc_fractions = YAXArray((classax,),[abundance[i] for i in classkeys])
    )
    
end