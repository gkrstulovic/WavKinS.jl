using WavKinS
using NCDatasets
using LaTeXStrings
using GLMakie



main_dir = "/home/vincent/Data/Wavkins/WKE_Stratified_Asymp/"
#main_dir = "/home/vlabarre/Bureau/data/Wavkins/WKE_Stratified_Asymp/"
data_dirs = readdir(main_dir)


for id in eachindex(data_dirs)

    if occursin("80", data_dirs[id]) && occursin("decayingGM", data_dirs[id])

    data_dir = main_dir * data_dirs[id] * "/"

    display(data_dir)

    ###################
    # Loading Dataset #
    ###################
    ds = NCDataset(data_dir * "WKE_Stratified_Asymp_data.nc", "r")
    dsGlobal = ds.group["Global"]
    dsSpectral = ds.group["Spectral"]
    
    
    ##########
    # Global #
    ##########
    t = dsGlobal["Times"][:]
    N = dsGlobal["N"][:]
    H = dsGlobal["H"][:]
    Disp = dsGlobal["Disp"][:]

    fig, ax = plot_1D_base!(t, N; xlabel=L"$t$", ylabel=L"$N$", title="Total wave action", xscale=identity, yscale=identity)
    save(data_dir * "Total_Wave_Action.png", fig)
    fig, ax = plot_1D_base!(t, H; xlabel=L"$t$", ylabel=L"$H$", title="Energy", xscale=identity, yscale=identity)
    save(data_dir * "Total_Wave_Energy.png", fig)
    fig, ax = plot_1D_base!(t, Disp; xlabel=L"$t$", ylabel=L"$\mathcal{D}$", title="Energy dissipation", xscale=identity, yscale=identity)
    save(data_dir * "Total_Energy_Dissipation.png", fig)
    

    ############
    # Spectrum #
    ############
    t = dsSpectral["TimesSP"][:]
    kh = dsSpectral["kh"]
    kz = dsSpectral["kz"]
    nk = dsSpectral["nk"]
    Ek = dsSpectral["Ek"]
    Pkh = dsSpectral["Pkh"]
    Pkz = dsSpectral["Pkz"]
    Nk = wave_spectrum_khkz(kh[1],kh[end],length(kh),kz[1],kz[end],length(kz))
    Run = Stratified_Asymp(Nk; interp_scheeme=WavKinS.bilin_interp_khkz, time_stepping_scheeme=WavKinS.RK2_step);
    its = 1:5:length(t)
    

    kk = [kh[1], kh[end]]
    fig, ax, cb = plot_2D_slices_khkz!(Run; ρ=Run.ω, dir=1, nslices=10, title=L"$e(k_h,k_z,t=0)$", zlims=(1e-3, 1e4))  
    record(fig, data_dir * "ek_slices_kh_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        plot_2D_slices_khkz!(Run; ρ=Run.ω, dir=1, nslices=10, fig=fig, ax=ax, cb=cb)
        ax.title = L"e(k_h,k_z,t=%$(round(t[it],digits=1)))"
        plot_theo!(ax, kk, 1e1*kk.^-2)
    end

    fig, ax, cb = plot_2D_slices_khkz!(Run; ρ=Run.ω, dir=2, nslices=10, title=L"$e(k_h,k_z,t=0)$", zlims=(1e-3, 1e4))  
    record(fig, data_dir * "ek_slices_kz_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        plot_2D_slices_khkz!(Run; ρ=Run.ω, dir=2, nslices=10, fig=fig, ax=ax, cb=cb)
        ax.title = L"e(k_h,k_z,t=%$(round(t[it],digits=1)))"
        plot_theo!(ax, kk, 1e1*kk.^-2)
    end

    
    fig, ax, hm = plot_heatmap_base!(kh,kz,Ek[:,:,1]; xlabel=L"k_h", ylabel=L"k_z", title=L"$e(k_h,k_z,t=0)$", zlims=(1e-5, 1e3))    
    record(fig, data_dir * "ek_vs_time.mp4", its; framerate = 24) do it
        plot_heatmap_base!(kh,kz,Ek[:,:,it]; fig, ax, hm)
        ax.title = L"e(k_h,k_z,t=%$(round(t[it],digits=1)))"
    end

    fig, ax, cb = plot_slices_ωkz_spectrum!(Run; ρ=Run.ω, nslices=20, zlabel="", title=L"$e(\omega_k,k_z,t=0)$", zlims=(1e-5, 1e3))    
    record(fig, data_dir * "ek_slices_omegakz_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        plot_slices_ωkz_spectrum!(Run; ρ=Run.ω, nslices=20, fig=fig, ax=ax, cb=cb)
        ax.title = L"e(\omega_k,k_z,t=%$(round(t[it],digits=1)))"
    end

    ωk = [1e-2, 1e1]
    fig, ax, cb = plot_slices_kzω_spectrum!(Run; ρ=Run.ω, nslices=20, zlabel="", title=L"$e(\omega_k,k_z,t=0)$", zlims=(1e-5, 1e3), ωlims=(1e-2, 1e1))  
    record(fig, data_dir * "ek_slices_kzomega_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        plot_slices_kzω_spectrum!(Run; ρ=Run.ω, nslices=20, fig=fig, ax=ax, cb=cb)
        ax.title = L"e(\omega_k,k_z,t=%$(round(t[it],digits=1)))"
        plot_theo!(ax, ωk, 1e1*ωk.^-2)
    end

    fig, ax, hm = plot_ωkz_spectrum!(Run; ρ=Run.ω, title=L"$e(\omega_k,k_z,t=0)$", zlims=(1e-5, 1e3))    
    record(fig, data_dir * "ek_omegakz_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        fig, ax, hm = plot_ωkz_spectrum!(Run; ρ=Run.ω, fig=fig, ax=ax, hm=hm)
        ax.title = L"e(\omega_k,k_z,t=%$(round(t[it],digits=1)))"
    end

    fig, ax, cb = plot_slices_ξkz_spectrum!(Run; ρ=Run.ω, nslices=20, zlabel="", title=L"$e(\xi_k,k_z,t=0)$", zlims=(1e-8, 1e3))    
    record(fig, data_dir * "ek_slices_xikz_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        plot_slices_ξkz_spectrum!(Run; ρ=Run.ω, nslices=20, fig=fig, ax=ax, cb=cb)
        ax.title = L"e(\xi_k,k_z,t=%$(round(t[it],digits=1)))"
    end

    ξk = [1e-1, 1e2]
    fig, ax, cb = plot_slices_kzξ_spectrum!(Run; ρ=Run.ω, nslices=20, zlabel="", title=L"$e(\xi_k,k_z,t=0)$", zlims=(1e-5, 1e3), ξlims=(1e-1, 1e2))  
    record(fig, data_dir * "ek_slices_kzxi_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        plot_slices_kzξ_spectrum!(Run; ρ=Run.ω, nslices=20, fig=fig, ax=ax, cb=cb)
        ax.title = L"e(\xi_k,k_z,t=%$(round(t[it],digits=1)))"
    end

    fig, ax, hm = plot_ξkz_spectrum!(Run; ρ=Run.ω, title=L"$e(\xi_k,k_z,t=0)$", zlims=(1e-5, 1e3))    
    record(fig, data_dir * "ek_xikz_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        fig, ax, hm = plot_ξkz_spectrum!(Run; ρ=Run.ω, fig=fig, ax=ax, hm=hm)
        ax.title = L"e(\xi_k,k_z,t=%$(round(t[it],digits=1)))"
    end

    fig, ax = plot_1D_spectra_density!(Run; ρ=Run.ω, ylabel="", title=L"$e(k_i,t=0)$", ylims=(1e-5, 1e3))    
    record(fig, data_dir * "1D_spectra_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        fig, ax = plot_1D_spectra_density!(Run; ρ=Run.ω, fig=fig, ax=ax)
        ax.title = L"e(k_i,t=%$(round(t[it],digits=1)))"
    end

    fig, ax = plot_ω_spectrum_density!(Run; ρ=Run.ω, ylabel="", title=L"$e(\omega_k,t=0)$", ωlims=(1e-2,1e1), ylims=(1e-5, 1e3))    
    record(fig, data_dir * "omega_spectrum_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        fig, ax = plot_ω_spectrum_density!(Run; ρ=Run.ω, fig=fig, ax=ax)
        ax.title = L"e(\omega_k,t=%$(round(t[it],digits=1)))"
    end


    
    ##########
    # Fluxes #
    ##########
    #=fig, ax = plot_2D_flux(Nk,Pkh[:,:,1],Pkz[:,:,1])   
    record(fig, data_dir * "fluxes_vs_time.mp4", its; framerate = 24) do it
        plot_2D_flux(Nk,Pkh[:,:,it],Pkz[:,:,it]) 
        ax.title = L"P_k(k_h,k_z,t=%$(round(t[it],digits=1)))"
    end=#
    
    #=fig, ax = plot_energy_flux_isotropic!(Run; ylims=(-0.001,0.001))    
    record(fig, data_dir * "fluxes_isotropic_vs_time.mp4", its; framerate = 24) do it
        @. Run.Nk.nk = nk[:,:,it]
        St_k!(Run)
        plot_energy_flux_isotropic!(Run; fig=fig, ax=ax)
    end=#
    @. Run.Nk.nk = nk[:,:,end]
    St_k!(Run)
    fig, ax = plot_energy_flux_isotropic!(Run)
    save(data_dir * "fluxes_isotropic.png", fig)    
    
    

    
    #######################
    # Energy conservation #
    #######################
    conservation_ratio = zeros(0)
    its = 1:100:length(t)

    for it in its
        @. Run.Nk.nk = nk[:,:,it]
        
        St_k!(Run)
        append!(conservation_ratio, energy_conservation_ratio(Run)) 
    end
    fig, ax = plot_1D_base!(t[its], conservation_ratio; xlabel=L"$t$", ylabel=L"$ECR$", title="Energy conservation ratio", xscale=identity, yscale=identity, ylims=(0,1))
    save(data_dir * "Energy_Conservation_Ratio.png", fig)
    

    end

end