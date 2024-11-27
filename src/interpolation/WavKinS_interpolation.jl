#############  linear interpolation #############
"""
    lin_interp_vs(k,k1,k2,n1,n2)

Return the linear interpolation of ``n_k`` at point `k` computed using the points `(k1,n1)` and `(k2,n2)`

``n_k = n_1 \\frac{k_2-k}{k_2-k_1} + n_2 \\frac{k-k_1}{k_2-k_1}.`` 

"""
function lin_interp_vs(k,k1,k2,n1,n2)
    return ( n1*(k2-k) + n2*(k-k1) ) / (k2 - k1)
end



@doc raw"""
    update_coeff_interp!(::lin_interp,Nk::wave_spectrum)

Compute the coefficients ``\alpha`` and ``\beta`` for linear interpolation such that ``n_{\bf k} = \alpha + \beta k``.
* `Nk`: Structure containing the wave action wave_spectrum
* `interp`: structure where the coefficients are stored

"""
function update_coeff_interp!(interp::lin_interp, Nk::wave_spectrum)
    nk = Nk.nk
    kk = Nk.kk
    λ = Nk.λ

    for i = 2:Nk.M
        interp.coeff.α[i] = -(nk[i] - λ * nk[i-1]) / (λ - 1.0)
        interp.coeff.β[i] = (nk[i] - nk[i-1]) / (kk[i-1] * (λ - 1.0))
    end
    if interp.coeff.α[2] < 0.0
        interp.coeff.α[1] = 0.0
        interp.coeff.β[1] = nk[1] / kk[1]
    else
        interp.coeff.α[1] = interp.coeff.α[2]
        interp.coeff.β[1] = interp.coeff.β[2]
    end
end


@doc raw"""
    val_nk(::lin_interp,Nk::waveaction,k)

Return the value of `Nk.nk` interpolated at `k` (with `lin_interp`). Coefficients should be updated.

See also [`update_coeff_interp!`](@ref)
"""
function val_nk(interp::lin_interp, Nk::wave_spectrum, k)
    if k > 0.0 && k <= Nk.kk[end]
        iInterp = ceil(Int, 1+log(k/Nk.kk[1]) / Nk.logλ)
        iInterp = clamp(iInterp, 1, Nk.M)
        nkInterp = interp.coeff.α[iInterp] + interp.coeff.β[iInterp] * k
    else
        nkInterp = 0.0
    end
    return max(nkInterp, 0.0)
end
########################################################



#############  linear interpolation in log scale #############

@doc raw"""
    update_coeff_interp!(::linlog_interp,Nk::wave_spectrum)

Compute the coefficients ``\alpha`` and ``C_0`` for linear interpolation in logarithmic scale (power law) such that ``n_{\bf k} = C_0 k^{-\alpha} = C_0 e^{- \alpha \log k}``.
* `Nk`: Structure containing the wave action wave_spectrum
* `interp`: structure where the coefficients are stored

"""
function update_coeff_interp!(interp::linlog_interp, Nk::wave_spectrum)
    nk = Nk.nk
    kk = Nk.kk
    logλ = Nk.logλ

    clean_waveaction!(Nk;min_nk=1e-100) # Needed to avoid log(0)
    for i = 2:Nk.M
        interp.coeff.α[i] = log(nk[i-1]/nk[i]) / logλ
        interp.coeff.C0[i] = nk[i] * exp(interp.coeff.α[i] * log(kk[i]))
    end
    interp.coeff.α[1] = interp.coeff.α[2]
    interp.coeff.C0[1] = interp.coeff.C0[2]
end


@doc raw"""
    val_nk(::linlog_interp,Nk::waveaction,k)

Return the value of `Nk.nk` interpolated at `k` (with `linlog_interp`). Coefficients should be updated.

See also [`update_coeff_interp!`](@ref)
"""
function val_nk(interp::linlog_interp, Nk::wave_spectrum, k)
    if k > Nk.kk[1] && k <= Nk.kk[end]
        iInterp = ceil(Int, 1+log(k/Nk.kk[1]) / Nk.logλ)
        iInterp = clamp(iInterp, 1, Nk.M)
        nkInterp = interp.coeff.C0[iInterp] * exp( - interp.coeff.α[iInterp] * log(k))
    else
        nkInterp = lin_interp_vs(k,Nk.kk[1],Nk.kk[2],Nk.nk[1],Nk.nk[2])
    end
    return max(nkInterp, 0.0)
end

########################################################



#############  power exponential interpolation #############

@doc raw"""
    update_coeff_interp!(::powexp_interp,Nk::wave_spectrum)

Compute the coefficients ``C_0``, ``\alpha`` and ``\beta`` for power law exponential interpolation such that ``n_{\bf k} = C_0 k^{-\alpha} e^{-\beta k}``.
* `Nk`: Structure containing the wave action wave_spectrum
* `interp`: structure where the coefficients are stored

"""
function update_coeff_interp!(interp::powexp_interp, Nk::wave_spectrum)
    M = Nk.M
    nk = Nk.nk
    kk = Nk.kk
    λ = Nk.λ
    logλ = Nk.logλ

    c0 = interp.coeff.c0
    β = interp.coeff.β
    α = interp.coeff.α

    Lognim = log(nk[1])
    Logni = log(nk[2])
    eps = 1.e-40
    for i = 2:(M-1)
        Lognip = log(nk[i+1])
        if nk[i-1] > eps && nk[i] > eps && nk[i+1] > eps
            α[i] = (Lognip + λ * Lognim - Logni * (1 + λ)) / ((λ - 1) * logλ)
            β[i] = -(Lognim - 2 * Logni + Lognip) / (kk[i-1] * (λ - 1)^2)
            c0[i] = nk[i] * kk[i]^(α[i]) * exp(β[i] * kk[i])
        else
            c0[i] = 0.0
            α[i] = 0.0
            β[i] = 0.0
        end

        Lognim = Logni
        Logni = Lognip

    end
    c0[1] = c0[2]
    α[1] = α[2]
    β[1] = β[2]
    c0[M] = c0[M-1]
    α[M] = α[M-1]
    β[M] = β[M-1]

end


@doc raw"""
    val_nk(::powexp_interp,Nk::waveaction,k)

Return the value of `Nk.nk` interpolated at `k` (with `powexp_interp`). Coefficients should be updated.

See also [`update_coeff_interp!`](@ref)
"""
function val_nk(interp::powexp_interp, Nk::wave_spectrum, k)
    c0 = interp.coeff.c0
    β = interp.coeff.β
    α = interp.coeff.α
    if k > 0.0 && k <= Nk.kk[end]
        iInterp = ceil(Int, 1+log(k/Nk.kk[1]) / Nk.logλ)
        iInterp = clamp(iInterp, 1, Nk.M)
        nkInterp = c0[iInterp] * k^(-α[iInterp]) * exp(-β[iInterp] * k)
    else
        nkInterp = 0.0
    end
    return max(nkInterp, 0.0)
end



########################################################

#############  power Gaussian interpolation #############

@doc raw"""
    update_coeff_interp!(::powexp_interp,Nk::wave_spectrum)

Compute the coefficients ``C_0``, ``\alpha`` and ``\beta`` for power law Gaussian interpolation such that ``n_{\bf k} = C_0 k^{-\alpha} e^{-\beta k^2}``.
* `Nk`: Structure containing the wave action wave_spectrum
* `interp`: structure where the coefficients are stored

"""
function update_coeff_interp!(interp::powGauss_interp, Nk::wave_spectrum)
    M = Nk.M
    nk = Nk.nk
    kk = Nk.kk
    λ = Nk.λ
    logλ = Nk.logλ

    c0 = interp.coeff.c0
    β = interp.coeff.β
    α = interp.coeff.α

    Lognim = log(nk[1])
    Logni = log(nk[2])
    eps = 1.e-40
    for i = 2:(M-1)
        Lognip = log(nk[i+1])
        if nk[i-1] > eps && nk[i] > eps && nk[i+1] > eps
            α[i] = (Lognip + λ^2 * Lognim - Logni * (1 + λ^2)) / ((λ^2 - 1) * logλ)
            β[i] = -(Lognim - 2 * Logni + Lognip) / (kk[i-1]^2 * (λ^2 - 1)^2)
            c0[i] = nk[i] * kk[i]^(α[i]) * exp(β[i] * kk[i]^2)
        else
            c0[i] = 0.0
            α[i] = 0.0
            β[i] = 0.0
        end

        Lognim = Logni
        Logni = Lognip

    end
    c0[1] = c0[2]
    α[1] = α[2]
    β[1] = β[2]
    c0[M] = c0[M-1]
    α[M] = α[M-1]
    β[M] = β[M-1]

end


@doc raw"""
    val_nk(interp::powGauss_interp, Nk::wave_spectrum, k)

Return the value of `Nk.nk` interpolated at `k` (with `powGauss_interp`). Coefficients should be updated.

See also [`update_coeff_interp!`](@ref)
"""
function val_nk(interp::powGauss_interp, Nk::wave_spectrum, k)
    c0 = interp.coeff.c0
    β = interp.coeff.β
    α = interp.coeff.α
    if k > 0.0 && k <= Nk.kk[end]
        iInterp = ceil(Int, 1+log(k/Nk.kk[1]) / Nk.logλ)
        iInterp = clamp(iInterp, 1, Nk.M)
        nkInterp = c0[iInterp] * k^(-α[iInterp]) * exp(-β[iInterp] * k^2)
    else
        nkInterp = 0.0
    end
    return max(nkInterp, 0.0)
end
########################################################


#############  BS spline interpolation #############

@doc raw"""
    update_coeff_interp!(interp::BS_interp,Nk::wave_spectrum)

Compute the coefficients for B-spline interpolation. The method is given by the [BSplineKit](https://jipolanco.github.io/BSplineKit.jl/dev/) package.
* `Nk`: Structure containing the wave action wave_spectrum
* `interp`: structure where the coefficients are stored
"""
function update_coeff_interp!(interp::BS_interp, Nk::wave_spectrum)
    interpolate!(interp.coeff.spline, Nk.nk)
end

function val_nk(interp::BS_interp, Nk::wave_spectrum, k)
    return interp.coeff.spline(k)
end

########################################################

#############  bilinear interpolation khkz #############
"""
    bilin_interp_vs(kh,kz,khinf,khsup,kzinf,kzsup,n1,n2,n3,n4)

Return the bilinear interpolation of ``n_{\\bf k}`` at point (`kh`,`kz`) computed using the points (`khinf`,`kzinf`,`n1`), (`khsup`,`kzinf`,`n2`), (`khsup`,`kzsup`,`n3`) and (`khinf`,`kzsup`,`n4`).

"""
function bilin_interp_vs(kh,kz,khinf,khsup,kzinf,kzsup,n1,n2,n3,n4)
    w1 = (khsup-kh)*(kzsup-kz)
    w4 = (khsup-kh)*(kz-kzinf)
    w2 = (kh-khinf)*(kzsup-kz)
    w3 = (kh-khinf)*(kz-kzinf)
    return ( w1*n1 + w2*n2 + w3*n3 + w4*n4 ) / ((khsup-khinf)*(kzsup-kzinf))
end

@doc raw"""
    update_coeff_interp!(::bilin_interp_khkz, Nk::wave_spectrum_khkz)

Compute the coefficients ``C_0``, ``\alpha_h``, ``\alpha_z`` and ``\beta`` for bilinear interpolation such that ``n_{\bf k} = C_0 - \alpha_h k_h - \alpha_z k_z - \beta k_h k_z``.
* `Nk`: Structure containing the wave action wave_spectrum
* `interp`: structure where the coefficients are stored

"""
function update_coeff_interp!(interp::bilin_interp_khkz, Nk::wave_spectrum_khkz)
    nk = Nk.nk
    kkh = Nk.kkh
    kkz = Nk.kkz
    λh = Nk.λh
    λz = Nk.λz

    for ih = 2:Nk.Mh, iz = 2:Nk.Mz
        n1 = nk[ih-1, iz-1]
        n2 = nk[ih, iz-1]
        n3 = nk[ih, iz]
        n4 = nk[ih-1, iz]
        b = (-n1 + n2 - n3 + n4) / ((λh - 1.0) * (λz - 1.0))

        interp.coeff.β[ih, iz] = b / (kkh[ih-1] * kkz[iz-1])
        interp.coeff.αh[ih, iz] = ((n1 - n2) / (λh - 1.0) - b) / kkh[ih-1]
        interp.coeff.αz[ih, iz] = ((n1 - n4) / (λz - 1.0) - b) / kkz[iz-1]
        interp.coeff.c0[ih, iz] = n1 + (n1 - n2) / (λh - 1.0) + (n1 - n4) / (λz - 1.0) - b
    end

    for iz = 2:Nk.Mz
        n1 = max(0.0, (kkh[2] * nk[1,iz-1] - kkh[1] * nk[2,iz-1]) / (kkh[2]-kkh[1]) )
        n2 = nk[1,iz-1]
        n3 = nk[1,iz]
        n4 = max(0.0, (kkh[2] * nk[1,iz] - kkh[1] * nk[2,iz]) / (kkh[2]-kkh[1]) )
        b = (-n1 + n2 - n3 + n4) / (kkh[1] * (kkz[iz] - kkz[iz-1]))
        
        interp.coeff.β[1, iz] = b
        interp.coeff.αh[1, iz] = (n1 - n2) / kkh[1] - b * kkz[iz-1] 
        interp.coeff.αz[1, iz] = (n1 - n4) / (kkz[iz] - kkz[iz-1])
        interp.coeff.c0[1, iz] = n1 + interp.coeff.αz[1, iz] * kkz[iz-1]
    end
    
    for ih = 2:Nk.Mh
        n1 = max(0.0, (kkz[2] * nk[ih-1,1] - kkz[1] * nk[ih-1,2]) / (kkz[2]-kkz[1]) )
        n2 = max(0.0, (kkz[2] * nk[ih,1] - kkz[1] * nk[ih,2]) / (kkz[2]-kkz[1]) )
        n3 = nk[ih,1]
        n4 = nk[ih-1,1]
        b = (-n1 + n2 - n3 + n4) / (kkz[1] * (kkh[ih] - kkh[ih-1]))

        interp.coeff.β[ih, 1] = b
        interp.coeff.αh[ih, 1] = (n1 - n2) / (kkh[ih] - kkh[ih-1])
        interp.coeff.αz[ih, 1] = (n1 - n4) / kkz[1] - b * kkh[ih-1] 
        interp.coeff.c0[ih, 1] = n1 + interp.coeff.αh[ih, 1] * kkh[ih-1]
    end

    # For ih = iz = 1 (linear nk = c0 - αh kh - αz kz if nks ≥ 0)
    n2 = max(0.0, (kkz[2] * nk[1,1] - kkz[1] * nk[1,2]) / (kkz[2]-kkz[1]) )
    n3 = nk[1,1]
    n4 = max(0.0, (kkh[2] * nk[1,1] - kkh[1] * nk[2,1]) / (kkh[2]-kkh[1]) )
    n1 = max(0.0, n2 - n3 + n4) 
    b = (-n1 + n2 - n3 + n4) / (kkz[1] * kkh[1]) 

    interp.coeff.β[1, 1] = b
    interp.coeff.αh[1, 1] = (n1 - n2) / kkh[1]
    interp.coeff.αz[1, 1] = (n1 - n4) / kkz[1]
    interp.coeff.c0[1, 1] = n1
end


@doc raw"""
    val_nk(::bilin_interp_khkz, Nk::wave_spectrum_khkz, kh, kz)

Return the value of `Nk.nk` interpolated at (`kh`, `kz`) (with `bilin_interp_khkz`). Coefficients should be updated.

See also [`update_coeff_interp!`](@ref)
"""
function val_nk(interp::bilin_interp_khkz, Nk::wave_spectrum_khkz, kh, kz)
    if kh > 0 && kh <= Nk.kkh[end] && kz > 0 && kz <= Nk.kkz[end]
        ihInterp = ceil(Int, 1+log(kh/Nk.kkh[1]) / Nk.logλh)
        ihInterp = clamp(ihInterp, 1, Nk.Mh)
        izInterp = ceil(Int, 1+log(kz/Nk.kkz[1]) / Nk.logλz)
        izInterp = clamp(izInterp, 1, Nk.Mz)
        nkInterp = interp.coeff.c0[ihInterp, izInterp] - interp.coeff.αh[ihInterp, izInterp] * kh - interp.coeff.αz[ihInterp, izInterp] * kz - interp.coeff.β[ihInterp, izInterp] * kh * kz
    else
        nkInterp = 0.0
    end
    return max(nkInterp, 0.0)
end
########################################################

#############  cpow interpolation khkz #############

@doc raw"""
    update_coeff_interp!(::cpow_interp_khkz, Nk::wave_spectrum_khkz)

Compute the coefficients ``C_0``, ``\alpha_h``, ``\alpha_z`` and ``\beta`` for bilinear interpolation such that ``n_{\bf k} = \beta + C_0  k_h^{-\alpha_h} k_z^{-\alpha_z k_z}``.
* `Nk`: Structure containing the wave action wave_spectrum
* `interp`: structure where the coefficients are stored

"""
function update_coeff_interp!(interp::cpow_interp_khkz, Nk::wave_spectrum_khkz)
    nk = Nk.nk
    kkh = Nk.kkh
    kkz = Nk.kkz
    logλh = Nk.logλh
    logλz = Nk.logλz

    for ih = 2:Nk.Mh, iz = 2:Nk.Mz
        n1 = nk[ih-1, iz-1]
        n2 = nk[ih, iz-1]
        n3 = nk[ih, iz]
        n4 = nk[ih-1, iz]
        b = (n1 * n3 - n2 * n4) / (n1 - n2 + n3 - n4)
        
        if b + 1e-20 < minimum([n1 n2 n3 n4]) 
            interp.coeff.β[ih, iz] = b
            interp.coeff.αh[ih, iz] = log((n1 - b) / (n2 - b)) / logλh
            interp.coeff.αz[ih, iz] = log((n1 - b) / (n4 - b)) / logλz
            interp.coeff.c0[ih, iz] = (n3 - b) * kkh[ih]^interp.coeff.αh[ih, iz] * kkz[iz]^interp.coeff.αz[ih, iz]
        else
            interp.coeff.β[ih, iz] = 0.25 * (n1 + n2 + n3 + n4)
            interp.coeff.αh[ih, iz] = 0
            interp.coeff.αz[ih, iz] = 0
            interp.coeff.c0[ih, iz] = 0
        end

    end

    # Not needed anymore because bilinear interpolation at the borders
    #=interp.coeff.β[1, :] = interp.coeff.β[2, :]
    interp.coeff.c0[1, :] = interp.coeff.c0[2, :]
    interp.coeff.αh[1, :] = interp.coeff.αh[2, :]
    interp.coeff.αz[1, :] = interp.coeff.αz[2, :]

    interp.coeff.β[:, 1] = interp.coeff.β[:, 2]
    interp.coeff.c0[:, 1] = interp.coeff.c0[:, 2]
    interp.coeff.αh[:, 1] = interp.coeff.αh[:, 2]
    interp.coeff.αz[:, 1] = interp.coeff.αz[:, 2]

    # Only power law interpolation for ih = iz = 1
    interp.coeff.β[1, 1] = 0.0
    interp.coeff.αh[1, 1] = interp.coeff.αh[2, 1]
    interp.coeff.αz[1, 1] = interp.coeff.αz[1, 2]
    interp.coeff.c0[1, 1] = nk[1, 1] * (kkh[1]^interp.coeff.αh[1, 1]) * (kkz[1]^interp.coeff.αz[1, 1])=# 
end


@doc raw"""
    val_nk(::cpow_interp_khkz, Nk::wave_spectrum_khkz, kh, kz)

Return the value of `Nk.nk` interpolated at (`kh`, `kz`) (with `cpow_interp_khkz`). Coefficients should be updated.

See also [`update_coeff_interp!`](@ref)
"""
function val_nk(interp::cpow_interp_khkz, Nk::wave_spectrum_khkz, kh, kz)
    if kh > 0 && kh <= Nk.kkh[end] && kz > 0 && kz <= Nk.kkz[end]
        ihInterp = ceil(Int, 1+log(kh/Nk.kkh[1]) / Nk.logλh)
        ihInterp = clamp(ihInterp, 1, Nk.Mh)
        izInterp = ceil(Int, 1+log(kz/Nk.kkz[1]) / Nk.logλz)
        izInterp = clamp(izInterp, 1, Nk.Mz)
        
        if ihInterp != 1 && izInterp != 1
            nkInterp = interp.coeff.β[ihInterp, izInterp] + interp.coeff.c0[ihInterp, izInterp] * (kh^(-interp.coeff.αh[ihInterp, izInterp])) * (kz^(-interp.coeff.αz[ihInterp, izInterp]))
        elseif ihInterp == 1 && izInterp != 1 # we evaluate n1 and n4 using linear interpolation first
            n2 = Nk.nk[1,izInterp-1]
            n3 = Nk.nk[1,izInterp]
            n1 = max(0.0,lin_interp_vs(0.0,Nk.kkh[1],Nk.kkh[2],Nk.nk[1,izInterp-1],Nk.nk[2,izInterp-1]))
            n4 = max(0.0,lin_interp_vs(0.0,Nk.kkh[1],Nk.kkh[2],Nk.nk[1,izInterp],Nk.nk[2,izInterp]))
            nkInterp = bilin_interp_vs(kh,kz,0.0,Nk.kkh[1],Nk.kkz[izInterp-1],Nk.kkz[izInterp],n1,n2,n3,n4)
        elseif ihInterp != 1 && izInterp == 1  # we evaluate n1 and n2 using linear interpolation first
            n3 = Nk.nk[ihInterp,1]
            n4 = Nk.nk[ihInterp-1,1]
            n1 = max(0.0,lin_interp_vs(0.0,Nk.kkz[1],Nk.kkz[2],Nk.nk[ihInterp-1,1],Nk.nk[ihInterp-1,2]))
            n2 = max(0.0,lin_interp_vs(0.0,Nk.kkz[1],Nk.kkz[2],Nk.nk[ihInterp,1],Nk.nk[ihInterp,2]))
            nkInterp = bilin_interp_vs(kh,kz,Nk.kkh[ihInterp-1],Nk.kkh[ihInterp],0.0,Nk.kkz[1],n1,n2,n3,n4)
        else # we evaluate n1, n2 and n4 using linear interpolation first
            n2 = max(0.0, (Nk.kkz[2] * Nk.nk[1,1] - Nk.kkz[1] * Nk.nk[1,2]) / (Nk.kkz[2]-Nk.kkz[1]) )
            n3 = Nk.nk[1,1]
            n4 = max(0.0, (Nk.kkh[2] * Nk.nk[1,1] - Nk.kkh[1] * Nk.nk[2,1]) / (Nk.kkh[2]-Nk.kkh[1]) )
            n1 = max(0.0, n2 - n3 + n4) 
            nkInterp = bilin_interp_vs(kh,kz,0.0,Nk.kkh[1],0.0,Nk.kkz[1],n1,n2,n3,n4)
        end

    else
        nkInterp = 0.0
    end
    return max(nkInterp, 0.0)
end
########################################################

#############  bilinear interpolation in log scale khkz #############

@doc raw"""
    update_coeff_interp!(::bilinlog_interp_khkz, Nk::wave_spectrum_khkz)

Compute the coefficients ``C_0``, ``\alpha_h``, ``\alpha_z`` and ``\beta`` for bilinear interpolation in logarithmic scale (power law) such that ``n_{\bf k} = e^{C_0 - \alpha_h \log k_h - \alpha_z \log k_z - \beta \log k_h \log k_z}``.
* `Nk`: Structure containing the wave action wave_spectrum
* `interp`: structure where the coefficients are stored

"""
function update_coeff_interp!(interp::bilinlog_interp_khkz, Nk::wave_spectrum_khkz)
    nk = Nk.nk
    kkh = Nk.kkh
    kkz = Nk.kkz
    logλh = Nk.logλh
    logλz = Nk.logλz

    clean_waveaction!(Nk;min_nk=1e-100) # Needed to avoid log(0)
    for ih = 2:Nk.Mh, iz = 2:Nk.Mz
        z1 = log(nk[ih-1, iz-1])
        z2 = log(nk[ih, iz-1])
        z3 = log(nk[ih, iz])
        z4 = log(nk[ih-1, iz])
        zz = -z1 + z2 - z3 + z4

        interp.coeff.β[ih, iz] = zz / (logλh * logλz)
        interp.coeff.αh[ih, iz] = (z1 - z2 - zz * (log(kkz[iz])/logλz - 1.0)) / logλh
        interp.coeff.αz[ih, iz] = (z1 - z4 - zz * (log(kkh[ih])/logλh - 1.0)) / logλz
        interp.coeff.c0[ih, iz] = z3 + interp.coeff.αh[ih, iz] * log(kkh[ih]) + interp.coeff.αz[ih, iz] * log(kkz[iz]) + interp.coeff.β[ih, iz] * log(kkh[ih]) * log(kkz[iz])
    end
    
    # Not needed anymore because bilinear interpolation at the borders
    #=interp.coeff.β[1, :] .= interp.coeff.β[2, :]
    interp.coeff.c0[1, :] .= interp.coeff.c0[2, :]
    interp.coeff.αh[1, :] .= interp.coeff.αh[2, :]
    interp.coeff.αz[1, :] .= interp.coeff.αz[2, :]

    interp.coeff.β[:, 1] .= interp.coeff.β[:, 2]
    interp.coeff.c0[:, 1] .= interp.coeff.c0[:, 2]
    interp.coeff.αh[:, 1] .= interp.coeff.αh[:, 2]
    interp.coeff.αz[:, 1] .= interp.coeff.αz[:, 2]

    interp.coeff.β[1, 1] = interp.coeff.β[2, 2]
    interp.coeff.c0[1, 1] = interp.coeff.c0[2, 2]
    interp.coeff.αh[1, 1] = interp.coeff.αh[2, 2]
    interp.coeff.αz[1, 1] = interp.coeff.αz[2, 2]=#

end

@doc raw"""
    val_nk(::bilinlog_interp_khkz, Nk::wave_spectrum_khkz, kh, kz)

Return the value of `Nk.nk` interpolated at (`kh`, `kz`) (with `bilinlog_interp_khkz`). Coefficients should be updated.

See also [`update_coeff_interp!`](@ref)
"""
function val_nk(interp::bilinlog_interp_khkz, Nk::wave_spectrum_khkz, kh, kz)
    if kh > 0 && kh <= Nk.kkh[end] && kz > 0 && kz <= Nk.kkz[end]
        ihInterp = ceil(Int, 1+log(kh/Nk.kkh[1]) / Nk.logλh)
        ihInterp = clamp(ihInterp, 1, Nk.Mh)
        izInterp = ceil(Int, 1+log(kz/Nk.kkz[1]) / Nk.logλz)
        izInterp = clamp(izInterp, 1, Nk.Mz)
        
        if ihInterp != 1 && izInterp != 1
            nkInterp = exp( interp.coeff.c0[ihInterp, izInterp] - interp.coeff.αh[ihInterp, izInterp] * log(kh) - interp.coeff.αz[ihInterp, izInterp] * log(kz) - interp.coeff.β[ihInterp, izInterp] * log(kh) * log(kz) )
        elseif ihInterp == 1 && izInterp != 1 # we evaluate n1 and n4 using linear interpolation first
            n2 = Nk.nk[1,izInterp-1]
            n3 = Nk.nk[1,izInterp]
            n1 = max(0.0,lin_interp_vs(0.0,Nk.kkh[1],Nk.kkh[2],Nk.nk[1,izInterp-1],Nk.nk[2,izInterp-1]))
            n4 = max(0.0,lin_interp_vs(0.0,Nk.kkh[1],Nk.kkh[2],Nk.nk[1,izInterp],Nk.nk[2,izInterp]))
            nkInterp = bilin_interp_vs(kh,kz,0.0,Nk.kkh[1],Nk.kkz[izInterp-1],Nk.kkz[izInterp],n1,n2,n3,n4)
        elseif ihInterp != 1 && izInterp == 1  # we evaluate n1 and n2 using linear interpolation first
            n3 = Nk.nk[ihInterp,1]
            n4 = Nk.nk[ihInterp-1,1]
            n1 = max(0.0,lin_interp_vs(0.0,Nk.kkz[1],Nk.kkz[2],Nk.nk[ihInterp-1,1],Nk.nk[ihInterp-1,2]))
            n2 = max(0.0,lin_interp_vs(0.0,Nk.kkz[1],Nk.kkz[2],Nk.nk[ihInterp,1],Nk.nk[ihInterp,2]))
            nkInterp = bilin_interp_vs(kh,kz,Nk.kkh[ihInterp-1],Nk.kkh[ihInterp],0.0,Nk.kkz[1],n1,n2,n3,n4)
        else # we evaluate n1, n2 and n4 using linear interpolation first
            n2 = max(0.0, (Nk.kkz[2] * Nk.nk[1,1] - Nk.kkz[1] * Nk.nk[1,2]) / (Nk.kkz[2]-Nk.kkz[1]) )
            n3 = Nk.nk[1,1]
            n4 = max(0.0, (Nk.kkh[2] * Nk.nk[1,1] - Nk.kkh[1] * Nk.nk[2,1]) / (Nk.kkh[2]-Nk.kkh[1]) )
            n1 = max(0.0, n2 - n3 + n4) 
            nkInterp = bilin_interp_vs(kh,kz,0.0,Nk.kkh[1],0.0,Nk.kkz[1],n1,n2,n3,n4)
        end
    else
        nkInterp = 0.0
    end
    return max(nkInterp, 0.0)
end
########################################################







#############  exponential interpolation khkz #############

@doc raw"""
    update_coeff_interp!(::exp_interp_khkz, Nk::wave_spectrum_khkz)

Compute the coefficients ``C_0``, ``\alpha_h``, ``\alpha_z`` and ``\beta`` for exponential interpolation such that ``n_{\bf k} = e^{C_0 - \alpha_h k_h - \alpha_z k_z - \beta k_h k_z}``.
* `Nk`: Structure containing the wave action wave_spectrum
* `interp`: structure where the coefficients are stored

"""
function update_coeff_interp!(interp::exp_interp_khkz, Nk::wave_spectrum_khkz)
    nk = Nk.nk
    kkh = Nk.kkh
    kkz = Nk.kkz
    λh = Nk.λh
    λz = Nk.λz

    clean_waveaction!(Nk;min_nk=1e-100) # Needed to avoid log(0)
    for ih = 2:Nk.Mh, iz = 2:Nk.Mz
        Δkh = kkh[ih] - kkh[ih-1]
        Δkz = kkz[iz] - kkz[iz-1]
        z1 = log(nk[ih-1, iz-1])
        z2 = log(nk[ih, iz-1])
        z3 = log(nk[ih, iz])
        z4 = log(nk[ih-1, iz])
        zz = -z1 + z2 - z3 + z4

        interp.coeff.β[ih, iz] = zz / (Δkh * Δkz)
        interp.coeff.αh[ih, iz] = (z1 - z2 - zz * kkz[iz] / (Δkz * λz) ) / Δkh
        interp.coeff.αz[ih, iz] = (z1 - z4 - zz * kkh[ih] / (Δkh * λh) ) / Δkz
        interp.coeff.c0[ih, iz] = z3 + interp.coeff.αh[ih, iz] * kkh[ih] + interp.coeff.αz[ih, iz] * kkz[iz] + interp.coeff.β[ih, iz] * kkh[ih] * kkz[iz]
    end
    
    # Not needed anymore because bilinear interpolation at the borders
    #=interp.coeff.β[1, :] .= interp.coeff.β[2, :]
    interp.coeff.c0[1, :] .= interp.coeff.c0[2, :]
    interp.coeff.αh[1, :] .= interp.coeff.αh[2, :]
    interp.coeff.αz[1, :] .= interp.coeff.αz[2, :]

    interp.coeff.β[:, 1] .= interp.coeff.β[:, 2]
    interp.coeff.c0[:, 1] .= interp.coeff.c0[:, 2]
    interp.coeff.αh[:, 1] .= interp.coeff.αh[:, 2]
    interp.coeff.αz[:, 1] .= interp.coeff.αz[:, 2]=#

end

@doc raw"""
    val_nk(::exp_interp_khkz, Nk::wave_spectrum_khkz, kh, kz)

Return the value of `Nk.nk` interpolated at (`kh`, `kz`) (with `exp_interp_khkz`). Coefficients should be updated.

See also [`update_coeff_interp!`](@ref)
"""
function val_nk(interp::exp_interp_khkz, Nk::wave_spectrum_khkz, kh, kz)
    if kh > 0 && kh <= Nk.kkh[end] && kz > 0 && kz <= Nk.kkz[end]
        ihInterp = ceil(Int, 1+log(kh/Nk.kkh[1]) / Nk.logλh)
        ihInterp = clamp(ihInterp, 1, Nk.Mh)
        izInterp = ceil(Int, 1+log(kz/Nk.kkz[1]) / Nk.logλz)
        izInterp = clamp(izInterp, 1, Nk.Mz)
        
        if ihInterp != 1 && izInterp != 1
            nkInterp = exp( interp.coeff.c0[ihInterp, izInterp] - interp.coeff.αh[ihInterp, izInterp] * kh - interp.coeff.αz[ihInterp, izInterp] * kz - interp.coeff.β[ihInterp, izInterp] * kh * kz )
        elseif ihInterp == 1 && izInterp != 1 # we evaluate n1 and n4 using linear interpolation first
            n2 = Nk.nk[1,izInterp-1]
            n3 = Nk.nk[1,izInterp]
            n1 = max(0.0,lin_interp_vs(0.0,Nk.kkh[1],Nk.kkh[2],Nk.nk[1,izInterp-1],Nk.nk[2,izInterp-1]))
            n4 = max(0.0,lin_interp_vs(0.0,Nk.kkh[1],Nk.kkh[2],Nk.nk[1,izInterp],Nk.nk[2,izInterp]))
            nkInterp = bilin_interp_vs(kh,kz,0.0,Nk.kkh[1],Nk.kkz[izInterp-1],Nk.kkz[izInterp],n1,n2,n3,n4)
        elseif ihInterp != 1 && izInterp == 1  # we evaluate n1 and n2 using linear interpolation first
            n3 = Nk.nk[ihInterp,1]
            n4 = Nk.nk[ihInterp-1,1]
            n1 = max(0.0,lin_interp_vs(0.0,Nk.kkz[1],Nk.kkz[2],Nk.nk[ihInterp-1,1],Nk.nk[ihInterp-1,2]))
            n2 = max(0.0,lin_interp_vs(0.0,Nk.kkz[1],Nk.kkz[2],Nk.nk[ihInterp,1],Nk.nk[ihInterp,2]))
            nkInterp = bilin_interp_vs(kh,kz,Nk.kkh[ihInterp-1],Nk.kkh[ihInterp],0.0,Nk.kkz[1],n1,n2,n3,n4)
        else # we evaluate n1, n2 and n4 using linear interpolation first
            n2 = max(0.0, (Nk.kkz[2] * Nk.nk[1,1] - Nk.kkz[1] * Nk.nk[1,2]) / (Nk.kkz[2]-Nk.kkz[1]) )
            n3 = Nk.nk[1,1]
            n4 = max(0.0, (Nk.kkh[2] * Nk.nk[1,1] - Nk.kkh[1] * Nk.nk[2,1]) / (Nk.kkh[2]-Nk.kkh[1]) )
            n1 = max(0.0, n2 - n3 + n4) 
            nkInterp = bilin_interp_vs(kh,kz,0.0,Nk.kkh[1],0.0,Nk.kkz[1],n1,n2,n3,n4)
        end
    else
        nkInterp = 0.0
    end
    return max(nkInterp, 0.0)
end
########################################################


