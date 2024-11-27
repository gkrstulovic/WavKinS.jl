# @doc raw"""
#     get_global_diagnostics!(Run::MMT)
    
# Compute default global diagnostics
    
# * `Run`: MMT WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
# This routine computes and store in `Run.diags.glob_diag` the default diagnostics:
# The current time, the total waveaction, the total energy, and the total dissipation.

# """
# function get_global_diagnostics!(Run::MMT)
#     get_global_diagnostics_default!(Run)
#     return nothing
# end

@doc raw"""
    store_spectral(Run::MMT)
    
Compute and store MMT spectral quantities
    
* `Run`: MMT WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.sp_store`  the waveaction, and the energy and waveaction flux spectra.

"""
function store_spectral!(Run::MMT)
    store_spectral_default!(Run)
    Run.diags.sp_store["Qk"] = hcat(Run.diags.sp_store["Qk"], zeros(Run.Nk.M))
    waveaction_flux!(@view(Run.diags.sp_store["Qk"][:, end]), Run)
    return nothing
end

@doc raw"""
    compute_spectral(Run)
    
Compute and store MMT current spectral quantities
    
* `Run`: MMT WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.sp_outs` the waveaction, energy spectra and their corresponding fluxes.

"""
function compute_spectral!(Run::MMT)
    compute_spectral_default!(Run)
    if Run.diags.sp_outs["Qk"].write_sp
        waveaction_flux!(Run.diags.sp_outs["Qk"].sp, Run)
    end
    return nothing
end
