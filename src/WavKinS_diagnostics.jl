@doc raw"""
    get_global_diagnostics(Run)
    
Compute default global diagnostics
    
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.glob_diag` the default diagnostics:
The current time, the total waveaction, the total energy, and the total dissipation.

This routine is typically called from `get_global_diagnostics(Run)` which might add additional diagnostics for each physical system.

"""
function get_global_diagnostics!(Run)
    get_global_diagnostics_default!(Run)
    return nothing
end

function get_global_diagnostics_default!(Run)
    AA = total_waveaction(Run)
    Energy = energy(Run)
    Disp = energy_dissipation(Run)
    push!(Run.diags.glob_diag["Times"].out, Run.t)
    push!(Run.diags.glob_diag["N"].out, AA)
    push!(Run.diags.glob_diag["H"].out, Energy)
    push!(Run.diags.glob_diag["Disp"].out, Disp)
    return nothing
end

@doc raw"""
    store_spectral!(Run)
    
Compute and store default spectral quantities
    
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.sp_store` the default spectral quantities: the waveaction and the energy fluc spectrum.
This routine is typically called from `store_spectral(Run)` which might add additional spectra for each physical system.

"""
function store_spectral!(Run)
    store_spectral_default!(Run)
    return nothing
end

function store_spectral_default!(Run)
    Run.diags.sp_store["nk"] = hcat(Run.diags.sp_store["nk"], Run.Nk.nk)
    Run.diags.sp_store["Pk"] = hcat(Run.diags.sp_store["Pk"], zeros(Run.Nk.M))
    energy_flux!(@view(Run.diags.sp_store["Pk"][:, end]), Run)
    return nothing
end

@doc raw"""
    compute_spectral!(Run)
    
Compute and store default current spectral quantities
    
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.sp_outs` the default spectral quantities: the waveaction, energy, and the energy fluc spectra.
This routine is typically called from `compute_spectral!(Run)` which might add additional spectra for each physical system.

"""
function compute_spectral!(Run)
    compute_spectral_default!(Run)
    return nothing
end
    
function compute_spectral_default!(Run)
    if Run.diags.sp_outs["nk"].write_sp
        Run.diags.sp_outs["nk"].sp .= Run.Nk.nk
    end
    if Run.diags.sp_outs["Ek"].write_sp
        energy_spectrum!(Run.diags.sp_outs["Ek"].sp, Run)
    end
    if Run.diags.sp_outs["Pk"].write_sp
        energy_flux!(Run.diags.sp_outs["Pk"].sp, Run)
    end
    return nothing
end

@doc raw"""
    total_waveaction(Run)

Compute total wave action ``\int n_{\bf k} \mathrm{d}{\bf k}``. 
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
"""
function total_waveaction(Run)
    return total_integral_density(Run)
end

@doc raw"""
    energy(Run)
    
Compute total energy ``\int \omega_{\bf k} n_{\bf k} \mathrm{d}{\bf k}``. 
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}`` and the dispersion relation ``\omega_{\bf k}``

"""
function energy(Run)
    return total_integral_density(Run, Run.ω)
end

@doc raw"""
    waveaction_dissipation(Run)
    
Compute total wave action dissipation of the system ``\int d_{\bf k} n_{\bf k} \mathrm{d}{\bf k}``.
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}`` and the dissipation coefficients ``d_{\bf k}``
"""
function waveaction_dissipation(Run)
    return density_dissipation(Run)
end

@doc raw"""
    energy_dissipation(Run)
    
Compute total energy dissipation of the system ``\int d_{\bf k} \omega_{\bf k} n_{\bf k} \mathrm{d}{\bf k}``.
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``, the dispersion relation ``\omega_{\bf k}`` and the dissipation coefficients ``d_{\bf k}``
"""
function energy_dissipation(Run)
    return density_dissipation(Run, Run.ω)
end

@doc raw"""
    energy_injection(Run)
    
Compute total energy injection of the system ``\int \omega_{\bf k} f_{\bf k} \mathrm{d}{\bf k}``.
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``, the dispersion relation ``\omega_{\bf k}`` and the forcing coefficients ``f_{\bf k}``
"""
function energy_injection(Run)
    return density_injection(Run, Run.ω)
end

@doc raw"""
    energy_spectrum!(Ek::Union{AbstractVector,AbstractMatrix},Run,ρ=one)

Compute energy spectrum. 
* `Ek`: where energy spectrum is stored
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
"""
function energy_spectrum!(Ek::Union{AbstractVector,AbstractMatrix}, Run)
    density_spectrum!(Ek, Run,Run.ω)
end

@doc raw"""
    waveaction_flux!(Qk::AbstractVector, Run)
    
Compute waveaction flux. For isotropic systems, it is ``P(k) = - \int\limits_{|{\bf k}'|<k} St_{{\bf k}'} \mathrm{d}{{\bf k}'}``.
* `Qk`: where energy flux is stored
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
"""
function waveaction_flux!(Qk::AbstractVector, Run)
    density_flux!(Qk, Run)
end

@doc raw"""
    energy_flux!(Pk::AbstractVector, Run)
    
Compute energy flux. For isotropic systems, it is ``P(k) = - \int\limits_{|{\bf k}'|<k} \omega_{{\bf k}'} St_{{\bf k}'} \mathrm{d}{{\bf k}'}``.
* `Pk`: where energy flux is stored
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}`` and dispersion relation ``\omega_{\rm k}``
"""
function energy_flux!(Pk::AbstractVector, Run)
    density_flux!(Pk, Run, Run.ω)
end

####### definition of general diagnostics for an invariant ρ  ##############@

function one(x...) 
    return 1
end


@doc raw"""
    total_integral_density(Run,ρ)
    
Compute total integral weighted by the density ``\rho`` as ``\int \rho_{\bf k} n_{\bf k} \mathrm{d}{\bf k}``.
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function total_integral_density(Run, ρ=one)
    d = Run.dimension
    dΩ = Run.dΩ

    if Run.Nk_arguments == 1
        kk = Run.Nk.kk
        for i in eachindex(kk)
            Run.F1.nk[i] = ρ(kk[i]) * Run.Nk.nk[i] * kk[i]^(d - 1) * dΩ
        end
        return integrate_with_log_bins(Run.F1)
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        for ih in eachindex(kkh), iz in eachindex(kkz)
            Run.F1.nk[ih, iz] = ρ(kkh[ih], kkz[iz]) * Run.Nk.nk[ih, iz] * kkh[ih]^(d - 1) * dΩ
        end
        return integrate(Run.integ_scheeme, Run.F1)
    else
        @error "density with three arguments not yet implemented"
    end
end

@doc raw"""
    density_dissipation(Run, ρ=one)
    
Compute total dissipation weighted by the density  ``\rho`` as  ``\int d_{\bf k} \rho_{\bf k} n_{\bf k} \mathrm{d}{\bf k}``.
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``and the dissipation coefficients ``d_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function density_dissipation(Run, ρ=one)
    d = Run.dimension
    dΩ = Run.dΩ
    D = Run.FD.D

    if Run.Nk_arguments == 1
        kk = Run.Nk.kk
        for i in eachindex(kk)
            Run.F1.nk[i] = D[i] * ρ(kk[i]) * Run.Nk.nk[i] * kk[i]^(d - 1) * dΩ
        end
        return integrate_with_log_bins(Run.F1)
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        for ih in eachindex(kkh), iz in eachindex(kkz)
            Run.F1.nk[ih, iz] = D[ih, iz] * ρ(kkh[ih], kkz[iz]) * Run.Nk.nk[ih, iz] * kkh[ih]^(d - 1) * dΩ
        end
        return integrate(Run.integ_scheeme, Run.F1)
    else
        @error "wave action with three arguments not yet implemented"
    end
end

@doc raw"""
    density_injection(Run, ρ=one)
    
Compute total injection weighted by the density  ``\rho`` as  ``\int \rho_{\bf k} f_{\bf k} \mathrm{d}{\bf k}``.    
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}`` and the forcing coefficients ``f_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function density_injection(Run, ρ=one)
    ω = Run.ω
    d = Run.dimension
    dΩ = Run.dΩ
    f = Run.FD.f

    if Run.Nk_arguments == 1
        kk = Run.Nk.kk
        for i in eachindex(kk)
            Run.F1.nk[i] = ω(kk[i]) * f[i] * kk[i]^(d - 1) * dΩ
        end
        return integrate_with_log_bins(Run.F1)
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        for ih in eachindex(kkh), iz in eachindex(kkz)
            Run.F1.nk[ih, iz] = ω(kkh[ih], kkz[iz]) * f[ih, iz] * kkh[ih]^(d - 1) * dΩ
        end
        return integrate(Run.integ_scheeme, Run.F1)
    else
        @error "wave action with three arguments not yet implemented"
    end
end

###### Definition of spectrum and fluxes
@doc raw"""
    density_spectrum!(Ek::Union{AbstractVector,AbstractMatrix}, Run, ρ=one)

Compute spectrum of a quantity with spectral density ``\rho``. For isotropic systems, it is ``s(k) = \rho_{\bf k} n_{\bf k} k^{d-1} \mathrm{d}\Omega``. For axisymmetric systems, it is ``s(k_h,k_z) = \rho_{\bf k} n_{\bf k} k_h^{d-1} \mathrm{d}\Omega``.
* `Ek`: where spectrum is stored
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``, physical dimension ``d`` and surface of the unit sphere (or of the isotropic sector) ``\mathrm{d}Ω``
* `ρ`: density. By default ρ()=1
"""
function density_spectrum!(Ek::Union{AbstractVector,AbstractMatrix},Run,ρ=one)
    d = Run.dimension
    dΩ = Run.dΩ

    if Run.Nk_arguments == 1
        kk = Run.Nk.kk
        for i in eachindex(kk)
            Ek[i] = ρ(kk[i]) * Run.Nk.nk[i] * kk[i]^(d - 1) * dΩ
        end
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        for ih in eachindex(kkh), iz in eachindex(kkz)
            Ek[ih, iz] = ρ(kkh[ih], kkz[iz]) * Run.Nk.nk[ih, iz] * kkh[ih]^(d - 1) * dΩ
        end
    else
        @error "wave action with three arguments not yet implemented"
    end
end

@doc raw"""
    isotropic_density_spectrum!(Run, ρ=one)

Compute isotropic spectrum of a quantity with spectral density ``\rho`` as ``s(k) = \frac{1}{\mathrm{d}k} \int\limits_{|{\bf k'}|=k}^{k+\mathrm{d}k} \rho_{\bf k} n_{\bf k} \mathrm{d}{{\bf k}'}``. 
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function isotropic_density_spectrum!(Run, ρ=one)

    if Run.Nk_arguments == 1
        kk = Run.Nk.kk
        Ek = zeros(Run.Nk.M)
        density_spectrum!(Ek,Run,ρ)
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        kmin = sqrt(kkh[1]^2 + kkz[1]^2)
        kmax = sqrt(kkh[end]^2 + kkz[end]^2)
        kk = LogRange(kmin,kmax,Run.Nk.Mh)
        Ekk = zeros(Run.Nk.Mh,Run.Nk.Mz)
        Ek = zeros(Run.Nk.Mh)
        f = zeros(Run.Nk.Mh)
        density_spectrum!(Ekk,Run,ρ)
        
        for ik in 2:length(kk)
            @. f = sqrt(max.(kk[ik]^2 - kkh.^2, 0.0));
            R = area_ratio_logbins(kkh,kkz,Run.Nk.logλz,f);
            @. f = sqrt(max.(kk[ik-1]^2 - kkh.^2, 0.0));
            R = R - area_ratio_logbins(kkh,kkz,Run.Nk.logλz,f);
            @. Run.F1.nk = Ekk .* R 
            Ek[ik] = integrate(Run.integ_scheeme, Run.F1) / (kk[ik]-kk[ik-1])
        end        
    else
        @error "wave action with three arguments not yet implemented"
    end

    return kk, Ek
end

@doc raw"""
    kh_density_spectrum!(Run, ρ=one)

For axisymmetric systems, compute ``k_h`` spectrum (after average over ``k_z``) of a quantity with spectral density ``\rho``. 
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function kh_density_spectrum!(Run, ρ=one)

    if Run.Nk_arguments == 1
        @error "No kh spectrum for isotropic systems. Use density_spectrum(Run; ρ)."
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        Ekk = zeros(Run.Nk.Mh,Run.Nk.Mz)
        Ek = zeros(Run.Nk.Mh)
        density_spectrum!(Ekk,Run,ρ)
        
        for ih in eachindex(kkh)
            Ek[ih] = integrate_with_log_bins(Run.Nk.Mz, Run.Nk.λz, Run.Nk.logλz, kkz, Ekk[ih,:])
        end        
    else
        @error "wave action with three arguments not yet implemented"
    end

    return Ek
end

@doc raw"""
    kz_density_spectrum!(Run, ρ=one)

For axisymmetric systems, compute ``k_z`` spectrum (after average over ``k_h``) of a quantity with spectral density ``\rho``. 
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function kz_density_spectrum!(Run, ρ=one)

    if Run.Nk_arguments == 1
        @error "No kz spectrum for isotropic systems. Use density_spectrum(Run; ρ)."
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        Ekk = zeros(Run.Nk.Mh,Run.Nk.Mz)
        Ek = zeros(Run.Nk.Mz)
        density_spectrum!(Ekk,Run,ρ)
        
        for iz in eachindex(kkz)
            Ek[iz] = integrate_with_log_bins(Run.Nk.Mh, Run.Nk.λh, Run.Nk.logλh, kkh, Ekk[:,iz])
        end        
    else
        @error "wave action with three arguments not yet implemented"
    end

    return Ek
end

@doc raw"""
     density_flux!(Pk::AbstractVector, Run,ρ=one)
    
Compute flux of a density ``\rho``. For isotropic systems, it is ``P(k) = - \int\limits_{|{\bf k}'|<k} \rho_{{\bf k}'} St_{{\bf k}'} \mathrm{d}{{\bf k}'}``.
* `Pk`: AbstractVector where flux is stored
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function density_flux!(Pk::AbstractVector, Run, ρ=one)
    d = Run.dimension
    dΩ = Run.dΩ

    if Run.Nk_arguments == 1
        kk = Run.Nk.kk
        for i in eachindex(kk)
            Run.F1.nk[i] = -ρ(kk[i]) * Run.Sk.nk[i] * kk[i]^(d - 1) * dΩ
        end
        cumintegrate!(Pk, Run.F1)
    elseif Run.Nk_arguments == 2
        Mh = Run.Nk.Mh
        Mz = Run.Nk.Mz
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        λh = Run.Nk.λh
        λz = Run.Nk.λz
        logλh = Run.Nk.logλh
        logλz = Run.Nk.logλz
        Pkh = Pk[1]
        Pkz = Pk[2]
        for ih in eachindex(kkh), iz in eachindex(kkz)
            Run.F1.nk[ih,iz] = -ρ(kkh[ih],kkz[iz]) * Run.Sk.nk[ih,iz] * kkh[ih]^(d - 1) * dΩ
        end

        # Computing sum for horizontal flux
        for iz in eachindex(kkz)
            cumintegrate!(@view(Pkh[:,iz]), Mh, λh, logλh, kkh, Run.F1.nk[:,iz])
        end
        
        # Computing sum for vertical flux
        for ih in eachindex(kkh)
            cumintegrate!(@view(Pkz[ih,:]), Mz, λz, logλz, kkz, Run.F1.nk[ih,:])
        end

        # Prefactor kh^(d - 1)
        for ih in eachindex(kkh)
            Pkh[ih,:] = Pkh[ih,:] / (2 * kkh[ih]^(d - 1))
            Pkz[ih,:] = Pkz[ih,:] / (2 * kkh[ih]^(d - 1))
        end

    else 
        @error "wave action with three arguments not yet implemented"
    end
end

@doc raw"""
    density_flux_isotropic!(Pk::AbstractVector, Run,ρ=one)
    
Compute isotric (after sum over angle) flux of a density ``\rho``. It is ``P(k) = - \int\limits_{|{\bf k}'|<k} \rho_{{\bf k}'} St_{{\bf k}'} \mathrm{d}{{\bf k}'}``.
* `Pk`: AbstractVector where flux is stored
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function density_flux_isotropic!(Pk::AbstractVector, Run, ρ=one)
    d = Run.dimension
    dΩ = Run.dΩ

    if Run.Nk_arguments == 1
        density_flux!(Pk, Run, ρ)
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        kk = Run.Nk.kkh
        d = Run.dimension
        dΩ = Run.dΩ
        f = zeros(Run.Nk.Mh)
    
        #TODO: Find smarter algorithm with cumsum? For now, it computes the complete sum for all kk... (not efficient) 
        for ik in eachindex(kk)
          @. f = sqrt(max.(kk[ik]^2 - kkh.^2, 0.0));
          R = area_ratio_logbins(kkh,kkz,Run.Nk.logλz,f);
          for ih in eachindex(kkh), iz in eachindex(kkz)
            Run.F1.nk[ih,iz] = -ρ(kkh[ih],kkz[iz]) * Run.Sk.nk[ih,iz] * kkh[ih]^(d - 1) * dΩ * R[ih,iz] 
          end
          Pk[ik] = integrate(Run.integ_scheeme, Run.F1)
        end
    else 
        @error "wave action with three arguments not yet implemented"
    end
end

@doc raw"""
    density_flux_angular!(thk::AbstractVector, Pthk::AbstractVector, Run, ρ=one)
    
Compute angular flux of a density ``\rho``. It is ``P(\theta_{\bf k}) = - \int\limits_{\theta_{{\bf k}'}<\theta_{{\bf k}}} \rho_{{\bf k}'} St_{{\bf k}'} \mathrm{d}{{\bf k}'}``.
* `thk`: AbstractVector for the angles grid
* `Pthk`: AbstractVector where flux is stored
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function density_flux_angular!(thk::AbstractVector, Pthk::AbstractVector, Run, ρ=one)
    kkh = Run.Nk.kkh
    kkz = Run.Nk.kkz
    d = Run.dimension
    dΩ = Run.dΩ

    f = zeros(Run.Nk.Mh)

    #TODO: Find smarter algorithm with cumsum? For now, it computes the complete sum for all kk... (not efficient) 
    for ik in eachindex(thk)
      @. f = kkh * tan(thk[ik]);
      R = area_ratio_logbins(kkh,kkz,Run.Nk.logλz,f);
      for ih in eachindex(kkh), iz in eachindex(kkz)
        Run.F1.nk[ih,iz] = -ρ(kkh[ih],kkz[iz]) * Run.Sk.nk[ih,iz] * kkh[ih]^(d - 1) * dΩ * R[ih,iz] 
      end
      Pthk[ik] = integrate(Run.integ_scheeme, Run.F1)
    end
end

@doc raw"""
    total_density_flux(Run,ρ)
    
Compute sum of the flux weighted by the density ``\rho`` as ``\int \rho_{\bf k} St_{\bf k} \mathrm{d}{\bf k}``.    
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}`` and the collision integral ``St_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
Note: Ideally, it should be zero for a dynamical invariant (e.g. energy should be conserved for Acoustic, Petviashvilli, Strat_Asymp, ...).
"""
function total_density_flux(Run, ρ=one)
    d = Run.dimension
    dΩ = Run.dΩ

    if Run.Nk_arguments == 1
        kk = Run.Nk.kk
        for i in eachindex(kk)
            Run.F1.nk[i] = ρ(kk[i]) * Run.Sk.nk[i] * kk[i]^(d - 1) * dΩ
        end
        return integrate_with_log_bins(Run.F1)
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        for ih in eachindex(kkh), iz in eachindex(kkz)
            Run.F1.nk[ih, iz] = ρ(kkh[ih], kkz[iz]) * Run.Sk.nk[ih, iz] * kkh[ih]^(d - 1) * dΩ
        end
        return integrate(Run.integ_scheeme, Run.F1)
    else
        @error "density with three arguments not yet implemented"
    end
end

@doc raw"""
    total_density_abs_flux(Run,ρ)
    
Compute sum of the absolute value of the flux weighted by the density ``\rho`` as ``\int \rho_{\bf k} |St_{\bf k}| \mathrm{d}{\bf k}``.
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}`` and the collision integral ``St_{\bf k}``
* `ρ`: a function defining the weight of the integral. By default ρ()=1
"""
function total_density_abs_flux(Run, ρ=one)
    d = Run.dimension
    dΩ = Run.dΩ

    if Run.Nk_arguments == 1
        kk = Run.Nk.kk
        for i in eachindex(kk)
            Run.F1.nk[i] = ρ(kk[i]) * abs(Run.Sk.nk[i]) * kk[i]^(d - 1) * dΩ
        end
        return integrate_with_log_bins(Run.F1)
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        for ih in eachindex(kkh), iz in eachindex(kkz)
            Run.F1.nk[ih, iz] = ρ(kkh[ih], kkz[iz]) * abs(Run.Sk.nk[ih, iz]) * kkh[ih]^(d - 1) * dΩ
        end
        return integrate(Run.integ_scheeme, Run.F1)
    else
        @error "density with three arguments not yet implemented"
    end
end

@doc raw"""
    energy_conservation_ratio(Run)
    
Compute the ratio ``\frac{\left| \int \omega_{\bf k} St_{\bf k} \mathrm{d}{\bf k} \right|}{\int \omega_{\bf k} |St_{\bf k}| \mathrm{d}{\bf k}}``.
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``, the dispersion relation ``\omega_{\bf k}`` and the collision integral ``St_{\bf k}``
Note: If the energy is a dynamical invariant, this ratio must tends to zero as the resolution increases (see [Eden *et al.*, J. Phys. Oceanogr. 49, 737-749 (2019)](https://journals.ametsoc.org/view/journals/phoc/49/3/jpo-d-18-0075.1.xml)).
"""
function energy_conservation_ratio(Run)
    d = Run.dimension
    dΩ = Run.dΩ

    total_energy_flux = total_density_flux(Run, Run.ω)
    total_energy_abs_flux = total_density_abs_flux(Run, Run.ω)

    return abs(total_energy_flux)/total_energy_abs_flux
end

@doc raw"""
    total_entropy(Run)

Compute total wave action entropy ``\int \log(n_{\bf k}) \mathrm{d}{\bf k}``.
* `Run`: WavKinS simulation structure containing the wave action ``n_{\bf k}``
"""
function total_entropy(Run)
    d = Run.dimension
    dΩ = Run.dΩ
    eps = 1e-80
    if Run.Nk_arguments == 1
        kk = Run.Nk.kk
        for i in eachindex(kk)
            Run.F1.nk[i] = log(Run.Nk.nk[i] + eps) * kk[i]^(d - 1) * dΩ
        end
        return integrate_with_log_bins(Run.F1)
    elseif Run.Nk_arguments == 2
        kkh = Run.Nk.kkh
        kkz = Run.Nk.kkz
        for ih in eachindex(kkh), iz in eachindex(kkz)
            Run.F1.nk[ih, iz] = log(Run.Nk.nk[ih, iz] + eps) * kkh[ih]^(d - 1) * dΩ
        end
        return integrate(Run.integ_scheeme, Run.F1)
    else
        @error "entropy with three arguments not yet implemented"
    end
end