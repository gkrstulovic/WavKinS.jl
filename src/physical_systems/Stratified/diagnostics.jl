@doc raw"""
    store_spectral(Run::Stratified_Asymp)
    
Compute and store Stratified_Asymp spectral quantities
    
* `Run`: Stratified_Asymp WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.sp_store`  the waveaction, and the energy and waveaction flux spectra.

"""
function store_spectral!(Run::Stratified_Asymp)
    Run.diags.sp_store["Pkh"] = cat(Run.diags.sp_store["Pkh"], zeros(Run.Nk.Mh,Run.Nk.Mz); dims=3)
    Run.diags.sp_store["Pkz"] = cat(Run.diags.sp_store["Pkz"], zeros(Run.Nk.Mh,Run.Nk.Mz); dims=3)
    Pkh = Run.diags.sp_store["Pkh"][:,:,end]
    Pkz = Run.diags.sp_store["Pkz"][:,:,end]
    energy_flux!([Pkh, Pkz], Run)
    return nothing
end

@doc raw"""
    compute_spectral(Run)
    
Compute and store Stratified_Asymp current spectral quantities
    
* `Run`: Stratified_Asymp WavKinS simulation structure containing the wave action ``n_{\bf k}``
    
This routine computes and store in `Run.diags.sp_outs` the waveaction, energy spectra and their corresponding fluxes.

"""
function compute_spectral!(Run::Stratified_Asymp)
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
