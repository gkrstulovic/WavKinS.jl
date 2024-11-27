@doc raw"""
    get_global_diagnostics!(Run::Petviashvilli_Asymp)
    
Compute default global diagnostics
    
* `Run`: Petviashvilli_Asymp simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.glob_diag` the default diagnostics:
The current time, the total waveaction, the total energy, total potential enstrophy, zonostrophy and the corresponding dissipations,

"""
function get_global_diagnostics!(Run::Petviashvilli_Asymp)
    get_global_diagnostics_default!(Run)
    Mx = total_integral_density(Run, ρ_Potential_Petviashvilli_Asymp)
    dMx = density_dissipation(Run, ρ_Potential_Petviashvilli_Asymp)
    Phi =  total_integral_density(Run,ρ_zonostrophy_Petviashvilli_Asymp)
    dPhi = density_dissipation(Run,ρ_zonostrophy_Petviashvilli_Asymp)

    push!(Run.diags.glob_diag["Mx"].out, Mx)
    push!(Run.diags.glob_diag["dMx"].out, dMx)
    push!(Run.diags.glob_diag["Phi"].out, Phi)
    push!(Run.diags.glob_diag["dPhi"].out, dPhi)
    return nothing
end

@doc raw"""
    store_spectral(Run::Union{Petviashvilli_Asymp,Petviashvilli})
    
Compute and store Petviashvilli_Asymp spectral quantities
    
* `Run`: Petviashvilli_Asymp or Petviashvilli WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.sp_store`  the waveaction, and the energy and waveaction flux spectra.

"""
function store_spectral!(Run::Union{Petviashvilli_Asymp,Petviashvilli})
    Run.diags.sp_store["Pkh"] = cat(Run.diags.sp_store["Pkh"], zeros(Run.Nk.Mh,Run.Nk.Mz); dims=3)
    Run.diags.sp_store["Pkz"] = cat(Run.diags.sp_store["Pkz"], zeros(Run.Nk.Mh,Run.Nk.Mz); dims=3)
    Pkh = Run.diags.sp_store["Pkh"][:,:,end]
    Pkz = Run.diags.sp_store["Pkz"][:,:,end]
    energy_flux!([Pkh, Pkz], Run)
    return nothing
end

@doc raw"""
    compute_spectral(Run)
    
Compute and store Petviashvilli_Asymp current spectral quantities
    
* `Run`: Petviashvilli_Asymp or Petviashvilli WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.sp_outs` the waveaction, energy spectra and their corresponding fluxes.

"""
function compute_spectral!(Run::Union{Petviashvilli_Asymp,Petviashvilli})
    if Run.diags.sp_outs["nk"].write_sp
        Run.diags.sp_outs["nk"].sp .= Run.Nk.nk
    end
    if Run.diags.sp_outs["Ek"].write_sp
        energy_spectrum!(Run.diags.sp_outs["Ek"].sp, Run)
    end
    if Run.diags.sp_outs["Pkh"].write_sp
        Pkh = Run.diags.sp_outs["Pkh"].sp
        Pkz = Run.diags.sp_outs["Pkz"].sp
        energy_flux!([Pkh,Pkz], Run)
    end
    return nothing
end
