@doc raw"""
    cumintegrate!(CumInt::AbstractVector, Sk::wave_spectrum) 

Compute cumulative integral of `Sk.Nk` using log bins and trapezoidal rule. The result is stored in `CumInt`.

"""
function cumintegrate!(CumInt::AbstractVector, Sk::wave_spectrum)
    M = Sk.M
    @assert length(CumInt) == M
    λ = Sk.λ
    logλ = Sk.logλ
    kk = Sk.kk
    F = Sk.nk

    cumintegrate!(CumInt, M, λ, logλ, kk, F)

end

@doc raw"""
    cumintegrate!(CumInt::AbstractVector, M::Int, λ::Float64, logλ::Float64, kk::Vector{Float64}, F::Vector{Float64}) 

Compute cumulative integral of `F` using log bins and trapezoidal rule. The result is stored in `CumInt`.
* `M`: number of grid points
* `λ`: logarithmic increment of the grid
* `logλ`: log(λ)
* `kk`: grid 

See also [`integrate_with_log_bins`](@ref).

"""
function cumintegrate!(CumInt::AbstractVector, M::Int, λ::Float64, logλ::Float64, kk::Vector{Float64}, F::Vector{Float64})
    @assert length(CumInt) == M

    CumInt .= 0.0
    CumInt[1] = integrate_with_log_bins(M, λ, logλ, kk, F, 1, 1)
    for i = 2:M
        CumInt[i] = CumInt[i-1] + integrate_with_log_bins(M, λ, logλ, kk, F, i, i)
    end

end

@doc raw"""
    cumintegraterev!(CumInt::AbstractVector, Sk::wave_spectrum) 

Compute cumulative integral using log bins in reverse order.

"""
function cumintegraterev!(CumInt::AbstractVector, Sk::wave_spectrum)
    M = Sk.M
    @assert length(CumInt) == M
    λ = Sk.λ
    logλ = Sk.logλ
    kk = Sk.kk
    F = Sk.nk

    cumintegraterevrev!(CumInt, M, λ, logλ, kk, F)

end

@doc raw"""
    cumintegraterev!(CumInt::AbstractVector, M::Int, λ::Float64, logλ::Float64, kk::Vector{Float64}, F::Vector{Float64}) 

Compute cumulative integral using log bins in reverse order.

"""
function cumintegraterev!(CumInt::AbstractVector, M::Int, λ::Float64, logλ::Float64, kk::Vector{Float64}, F::Vector{Float64})
    @assert length(CumInt) == M

    CumInt .= 0.0
    CumInt[M] = integrate_with_log_bins(M, λ, logλ, kk, F, M, M)
    for i = 1:M-1
        CumInt[M-i] = CumInt[M-i+1] + integrate_with_log_bins(M, λ, logλ, kk, F, M - i, M - i)
    end

end

@doc raw"""
    integrate_with_log_bins(Nk::wave_spectrum, imin::Int=1, imax::Int=-1)

Return the integral of `Nk.nk` using log bins and trapezoidal rule from index `imin` to `imax`.
By default, it integrates over the whole domain.

"""
function integrate_with_log_bins(Nk::wave_spectrum, imin::Int=1, imax::Int=-1)
    kk = Nk.kk
    M = Nk.M
    F = Nk.nk
    λ = Nk.λ
    logλ = Nk.logλ

    integrate_with_log_bins(M, λ, logλ, kk, F, imin, imax)
end

@doc raw"""
    integrate_with_log_bins(M::Int, λ::Float64, logλ::Float64, kk::Vector{Float64}, F::Vector{Float64}, imin::Int=1, imax::Int=-1)   

Return the integral of `F` using log bins and trapezoidal rule from index `imin` to `imax`.
* `M`: number of grid points
* `λ`: logarithmic increment of the grid
* `logλ`: log(λ)
* `kk`: grid 

By default, it integrates over the whole domain.
Namely ``I = \int\limits_{k_{\rm min}}^{k_{\rm max}} ~ F(k) ~ \mathrm{d}k = \int\limits_{i_{\rm min}}^{i_{\rm max}} ~ F(k=\lambda^{i} ) \lambda^{i} \ln \lambda ~ \mathrm{d} i``. It uses trapezoidal rule.

"""
function integrate_with_log_bins(M::Int, λ::Float64, logλ::Float64, kk::Vector{Float64}, F::Vector{Float64}, imin::Int=1, imax::Int=-1)
    F0 = 0.0
    Inte = 0.0

    if imax == -1
        imax = M
    end

    for i = imin:min(M, imax)
        if i == 1  # Fist domain is different
            F0 = max(0.0, -(F[2] - λ * F[1]) / (λ - 1.0)) #Linear extrapolation to k=0
            Inte += 0.5 * kk[1] * (F[1] + F0)
        else
            Inte += 0.5 * (F[i] * kk[i] + F[i-1] * kk[i-1]) * logλ
        end
    end

    return Inte

end

@doc raw"""
    integrate_with_log_bins_segment(Fa, Fb, ka, kb)

Return the integral of `F` using trapezoidal log bins from `ka` to `kb`. 
* `Fa`: value of `F` at `ka`
* `Fb`: value of `F` at `kb`

"""
function integrate_with_log_bins_segment(Fa, Fb, ka, kb)
    if ka != 0
        Inte = 0.5 * (Fa * ka + Fb * kb) * log(kb / ka)
    else
        Inte = 0.5 * (Fa + Fb) * kb
    end
    return Inte

end

@doc raw"""
    integrate_quad_with_log_bins(Nk::wave_spectrum, imin::Int=1, imax::Int=-1)

Return the integral of `Nk.nk` using log bins and Simpson rule from index `imin` to `imax`. 
By default, it integrates over the whole domain.

"""
function integrate_quad_with_log_bins(Sk::wave_spectrum, imin::Int=1, imax::Int=-1)
    kk = Sk.kk
    M = Sk.M
    F = Sk.nk
    logλ = Sk.logλ
    λ = Sk.λ
    F0 = 0.0
    Inte = 0.0

    if imax == -1
        imax = M
    end

    for i = imin:min(M, imax)
        if i == 1  # Fist domain is different
            F0 = max(0.0, -(F[2] - λ * F[1]) / (λ - 1.0)) #Linear extrapolation to k=0
            Inte += 0.5 * kk[1] * (F[1] + F0)
        elseif i < M
            Inte += 5 * F[i-1] * kk[i-1] * logλ / 12.0 + 2 * F[i] * kk[i] * logλ / 3.0 - F[i+1] * kk[i+1] * logλ / 12.0
        else
            # Inte += 5 * F[i-1]*kk[i-1] *logλ /12. +  2 * F[i]*kk[i] *logλ /3.
            Inte += 0.5 * (F[i] * kk[i] + F[i-1] * kk[i-1]) * logλ
        end
    end

    return Inte

end

@doc raw"""
    integrate(Nk::wave_spectrum, imin::Int=1, imax::Int=-1)

Return the integral of `Nk.nk` using (simple) trapezoidal rule from index `imin` to `imax`. 
By default, it integrates over the whole domain.

"""
function integrate(Sk::wave_spectrum, imin::Int=1, imax::Int=-1)
    kk = Sk.kk
    M = Sk.M
    F = Sk.nk
    λ = Sk.λ

    Inte = 0.0
    if imax == -1
        imax = M
    end

    for i = imin:min(M, imax)
        if i == 1  # Fist domain is different
            F0 = max(0.0, -(F[2] - λ * F[1]) / (λ - 1.0)) #Linear extrapolation to k=0
            Inte += 0.5 * kk[1] * (F[1] + F0)
        else
            Inte += 0.5 * (kk[i] - kk[i-1]) * (F[i] + F[i-1])
        end
    end

    return Inte

end

@doc raw"""
    integrate_with_grid(x, F, imin::Int=1, imax::Int=-1)

Return the integral of `F` using (simple) trapezoidal rule on grid `x` between index `imin` and `imax`. 

"""
function integrate_with_grid(x, F, imin::Int=2, imax::Int=-1)
    #@assert length(x) == length(F)
    Inte = 0.0

    if imax == -1
        imax = length(F)
    end

    for i in imin:min(length(x), imax)
        Inte += 0.5 * (x[i] - x[i-1]) * (F[i] + F[i-1])
    end
    return Inte
end

@doc raw"""
    integrate_with_grid(Sk::field_grid_1D, imin::Int=1, imax::Int=-1)

Return the integral of `Sk.F` using (simple) trapezoidal rule on grid `Sk.kk` between index `imin` and `imax`. 

"""
function integrate(Sk::field_grid_1D, imin::Int=2, imax::Int=-1)
    kk = Sk.kk
    M = Sk.M
    F = Sk.F

    return integrate_with_grid(kk, F, imin, imax)

end

"""
    clean_waveaction!(Nk::waveaction)

Set the minimum of the waveaction spectrum `Nk.nk` to `min_nk`. By default `min_nk=0` so it removes negative values.

"""
function clean_waveaction!(Nk::Union{wave_spectrum,wave_spectrum_khkz}; min_nk=0.0)

    for i in eachindex(Nk.nk)
        if Nk.nk[i] <= min_nk
            Nk.nk[i] = min_nk
        end
    end
end

@doc raw"""
    integrate(integ::integrate_with_log_bins_khkz, Nk::wave_spectrum_khkz, ihmin::Int=1, ihmax::Int=-1, izmin::Int=1, izmax::Int=-1)

Return the integral of `Nk.nk` using log bins and trapezoidal rule from index `ihmin` to `ihmax` and `izmin` to `izmax`. 
By default, it integrates over the whole domain.

Namely ``I = \int\limits_{k_{h{\rm min}}}^{k_{h{\rm max}}} \int\limits_{k_{z{\rm min}}}^{k_{z{\rm max}}} ~ F(k_h,k_z) ~ \mathrm{d}k_h \mathrm{d}k_z = \int\limits_{i_{h{\rm min}}}^{i_{h{\rm max}}} \int\limits_{i_{z{\rm min}}}^{i_{z{\rm max}}} ~ F(k_h=\lambda_h^{i_h}, k_z=\lambda_z^{i_z}) ~ \lambda_h^{i_h} \lambda_z^{i_z} \ln \lambda_h \ln \lambda_z ~ \mathrm{d} i_h \mathrm{d} i_z``. 

"""
function integrate(integ::integrate_with_log_bins_khkz, Nk::wave_spectrum_khkz, ihmin::Int=1, ihmax::Int=-1, izmin::Int=1, izmax::Int=-1)
    kkh = Nk.kkh
    kkz = Nk.kkz
    Mh = Nk.Mh
    Mz = Nk.Mz
    F = Nk.nk
    logλh = Nk.logλh
    logλz = Nk.logλz
    λh = Nk.λh
    λz = Nk.λz
    Inte = 0.0

    if ihmax == -1
        ihmax = Mh
    end
    if izmax == -1
        izmax = Mz
    end

    for ih = ihmin:min(Mh, ihmax), iz = izmin:min(Mz, izmax)
        if ih == 1 && iz == 1  # Fist domain is different
            F2 = max(0.0, F[1, 1] - (F[1, 2] - F[1, 1]) / (λz - 1.0))
            F4 = max(0.0, F[1, 1] - (F[2, 1] - F[1, 1]) / (λh - 1.0))
            Inte += 0.25 * (F2 + F4) * (kkh[1] * kkz[1]) #Linear extrapolation to k=0
        elseif ih == 1
            F1 = max(0.0, F[1, iz-1] - (F[2, iz-1] - F[1, iz-1]) / (λh - 1.0))
            F2 = F[1, iz-1]
            F3 = F[1, iz]
            F4 = max(0.0, F[1, iz] - (F[2, iz] - F[1, iz]) / (λh - 1.0))
            Inte += 0.25 * ((F1 + F2) * kkz[iz-1] + (F3 + F4) * kkz[iz]) * logλz * kkh[1] #Linear extrapolation in kh to kh=0
        elseif iz == 1
            F1 = max(0.0, F[ih-1, 1] - (F[ih-1, 2] - F[ih-1, 1]) / (λz - 1.0))
            F2 = max(0.0, F[ih, 1] - (F[ih, 2] - F[ih, 1]) / (λz - 1.0))
            F3 = F[ih, 1]
            F4 = F[ih-1, 1]
            Inte += 0.25 * ((F1 + F4) * kkh[ih-1] + (F2 + F3) * kkh[ih]) * logλh * kkz[1] #Linear extrapolation in kh to kz=0
        else
            Inte += 0.25 * (F[ih-1, iz-1] * kkh[ih-1] * kkz[iz-1] + F[ih, iz-1] * kkh[ih] * kkz[iz-1] + F[ih, iz] * kkh[ih] * kkz[iz] + F[ih-1, iz] * kkh[ih-1] * kkz[iz]) * logλh * logλz #trapezoidal integration with log bins
        end
    end

    return Inte

end

@doc raw"""
    integrate(integ::integrate_with_cpow_khkz, Nk::wave_spectrum_khkz, ihmin::Int=1, ihmax::Int=-1, izmin::Int=1, izmax::Int=-1)

Return the integral of `Nk.nk` assuming ``n_{\bf k} = \beta + C_0  k_h^{-\alpha_h} k_z^{-\alpha_z k_z}`` in each mesh from index `ihmin` to `ihmax` and `izmin` to `izmax`. 
By default, it integrates over the whole domain.

See also `cpow_interp_khkz`.

"""
function integrate(integ::integrate_with_cpow_khkz, Nk::wave_spectrum_khkz, ihmin::Int=1, ihmax::Int=-1, izmin::Int=1, izmax::Int=-1)
    kkh = Nk.kkh
    kkz = Nk.kkz
    Mh = Nk.Mh
    Mz = Nk.Mz
    F = Nk.nk
    logλh = Nk.logλh
    logλz = Nk.logλz
    λh = Nk.λh
    λz = Nk.λz
    Inte = 0.0

    interp = WavKinS.cpow_interp_khkz(kkh, kkz)
    WavKinS.update_coeff_interp!(interp, Nk)

    if ihmax == -1
        ihmax = Mh
    end
    if izmax == -1
        izmax = Mz
    end

    for ih = ihmin:min(Mh, ihmax), iz = izmin:min(Mz, izmax)
        if ih == 1 && iz == 1  # Fist domain is different
            F2 = max(0.0, F[1, 1] - (F[1, 2] - F[1, 1]) / (λz - 1.0))
            F4 = max(0.0, F[1, 1] - (F[2, 1] - F[1, 1]) / (λh - 1.0))
            Inte += 0.5 * (F2 + F4) * (kkh[1] * kkz[1]) #Linear extrapolation to k=0
        elseif ih == 1
            F1 = max(0.0, F[1, iz-1] - (F[2, iz-1] - F[1, iz-1]) / (λh - 1.0))
            F2 = F[1, iz-1]
            F3 = F[1, iz]
            F4 = max(0.0, F[1, iz] - (F[2, iz] - F[1, iz]) / (λh - 1.0))
            Inte += 0.25 * ((F1 + F2) * kkz[iz-1] + (F3 + F4) * kkz[iz]) * logλz * kkh[1] #Linear extrapolation in kh to kh=0
        elseif iz == 1
            F1 = max(0.0, F[ih-1, 1] - (F[ih-1, 2] - F[ih-1, 1]) / (λz - 1.0))
            F2 = max(0.0, F[ih, 1] - (F[ih, 2] - F[ih, 1]) / (λz - 1.0))
            F3 = F[ih, 1]
            F4 = F[ih-1, 1]
            Inte += 0.25 * ((F1 + F4) * kkh[ih-1] + (F2 + F3) * kkh[ih]) * logλh * kkz[1] #Linear extrapolation in kh to kz=0
        else
            c0 = interp.coeff.c0[ih, iz]
            αh = interp.coeff.αh[ih, iz]
            αz = interp.coeff.αz[ih, iz]
            β = interp.coeff.β[ih, iz]
            a = β * (1.0 - 1.0 / λh) * (1.0 - 1.0 / λz) * kkh[ih] * kkz[iz] + c0 * ((1.0 - λh^(αh - 1.0)) * (kkh[ih]^(1.0 - αh)) / (1.0 - αh)) * ((1.0 - λz^(αz - 1.0)) * (kkz[iz]^(1.0 - αz)) / (1.0 - αz))
            Inte += a
        end
    end

    return Inte

end