push!(LOAD_PATH, "../")
using WavKinS
using GLMakie; GLMakie.activate!()
using TimerOutputs
using Trapz

## Init basics
if Threads.nthreads() == 1
    @warn "Set the number of threads by doing export JULIA_NUM_THREADS "
end

const to = TimerOutput()

function run_simple_evol(Run, Param, Tfinal)
    step_scheeme = Run.time_stepping
    init_temporal_scheeme!(step_scheeme, Run, Param.dt)
    tini = Run.t
    tplotEvol = tini
    tglobalEvol = tini
    tsptoreEvol = tini
    tspwriteEvol = tini
    t = tini

    β = Run.β
    α = 1.0 / 2.0
    xQ = 2 * β / 3 + 1 - α / 3
    xP = 2 * β / 3 + 1
    #IpxP= dxI(xP,Q1,dIvec; β=0.)
    IpxP = WavKinS.dxI_MMT(xP; β)
    IpxQ = WavKinS.dxI_MMT(xQ; β)
    iQ = findfirst(x -> x >= 4 * 10^(-5), Run.Nk.kk)
    iP = findfirst(x -> x >= 0.015, Run.Nk.kk)
    PrefactorQ = 0.125
    PrefactorP = 0.125
    xcomp = 0


    @time while t < tini + Tfinal
        if any(isnan.(Run.Nk.nk))
            @warn "We got NaNs. Stop the temporal loop"
            break
        end
        @timeit to "Advance" advance!(step_scheeme, Run, Param.dt)
        t = Run.t

        @timeit to "Global quantities" if t >= tglobalEvol
            get_global_diagnostics!(Run)
            if Param.write_global
                output_global(Run, Param.outputDir)
            end
            Energy = Run.diags.glob_diag["H"].out[end]
            WaveAction = Run.diags.glob_diag["N"].out[end]
            println("t =", Run.t, "   H(t) =", Energy, "     N(t) = ", WaveAction)
            tglobalEvol += Param.tglobal
        end

        @timeit to "Plotting" if t >= tplotEvol
            plot_wave_action!(Run; fig=figN, ax=axN)
            plot_wave_action_flux!(Run; fig=figP, ax=axP)
            if size(Run.diags.sp_store["Pk"], 2) >= 1
                P0 = Run.diags.sp_store["Pk"][iP, end]
                Q0 = abs(Run.diags.sp_store["Qk"][iQ, end])
                PrefactorP = 2*abs(3.0 * P0 / (8 * pi * IpxP))^(1.0 / 3.0)
                PrefactorQ = 2*abs(3.0 * Q0 / (8 * pi * IpxQ))^(1.0 / 3.0)
                plot_theo!(axN, kk, PrefactorP * kk .^ (-xP + xcomp))
                plot_theo!(axN, kk, PrefactorQ * kk .^ (-xQ + xcomp); color="green")
            end
            sleep(0.0001)
            tplotEvol += Param.tplot
        end

        @timeit to "Storing spectral quantities" if t >= tsptoreEvol
            store_spectral!(Run)
            tsptoreEvol += Param.tspstore
        end

        @timeit to "Writing spectral quantities" if t >= tspwriteEvol
            if Param.write_spectral
                compute_spectral!(Run)
                output_spectra(Run, Param.outputDir)
            end
            tspwriteEvol += Param.tspwrite
        end


    end

end


###########################################################################################################
###########################################################################################################
## Defining the simulation

M = 256
kmin = 1e-6
kmax = 1e0

Nk = wave_spectrum(kmin,kmax,M)

# The MMT structure contains all the methods for integration, see help ?MMT for options
Run = MMT(Nk; interp_scheeme=WavKinS.lin_interp, drive_scheeme=WavKinS.RK2_step, β=0.0);



###########################################################################################################
##########################            running parmeters                           #########################
###########################################################################################################

isrestart = false # set this to true if simulations uses a restart

Tfinal = 8; # Final time of the simulation
dt = 0.001; # time step of the simulation

tplot = 0.01; # Plot every tplot times
tglobal = 0.01; # Compute, store and write global quantities every tglobal times 
tspstore = 0.01 # Compute and store spectral quantities every tspstore times

tspwrite = 0.1# write spectra every tspwrite
outputDir = "./" #output directory
write_global = true #write global
write_nk = true #write wave action
write_Ek = true #write energy spectrum
write_Pk = true #write energy flux  
write_Qk = true #write wave action flux   


# Creating simulations parameters and setting write_sp for each spectrum
write_spectral = write_nk || write_Ek || write_Pk || write_Qk
Param = simulation_parameters(dt, tplot, tglobal, tspstore, tspwrite,
    outputDir, write_global, write_spectral)

Run.diags.sp_outs["nk"].write_sp = write_nk
Run.diags.sp_outs["Ek"].write_sp = write_Ek
Run.diags.sp_outs["Pk"].write_sp = write_Pk
Run.diags.sp_outs["Qk"].write_sp = write_Qk


###########################################################################################################
###############     Setting forcing, dissipation and initial condition            #########################
###########################################################################################################


kk = Nk.kk
nu = 1
kd = 0.85 * kmax
kdLS = 3 * kmin
lapPower = 2
LSInvlapPower = 3
kf = 6.e-4;
Δkf = 5e-5


@. Run.Nk.nk = kk^2 * exp(-0.5 * (kk - kf)^2 / (Δkf)^2);
AA = total_waveaction(Run)
@. Run.Nk.nk = Run.Nk.nk / AA
@. Run.FD.f = Run.Nk.nk
@. Run.FD.D = nu * ((kk / kd)^(2 * lapPower) + (kk / kdLS)^(-LSInvlapPower))


if ~isrestart
    get_global_diagnostics!(Run)
    init_IO(Run, Param.outputDir)
    if Param.write_global
        output_global(Run, Param.outputDir)
    end
    if Param.write_spectral
        compute_spectral!(Run)
        output_spectra(Run, Param.outputDir)
    end
else
    load_spectrum(Run, Param.outputDir)
    get_global_diagnostics!(Run)
    compute_spectral!(Run)
end


#@. Run.Nk.nk = Run.Nk.nk /Eini
println("Run in resolultion M=", M, " with λ=", Run.Nk.λ)
println("Total wave action =", Run.diags.glob_diag["N"].out[end])
println("Total energy =", Run.diags.glob_diag["H"].out[end])


figN, axN = plot_wave_action!(Run; ylims=(1e-2, 1e+5))
figQ, axQ = plot_wave_action_flux!(Run; ylims=(-5, 5))
figP, axP = plot_energy_flux!(Run; ylims=(-5, 5))
display(figN) # figN (figP) if you want to plot wave action spectrum (wave action flux)
# #scatter!(ax, kk, Nk.nk)



## ######################################################################################################
##########################            Running the code                          #########################
#########################################################################################################

run_simple_evol(Run, Param, 0.0) # Good to run first to measure time properly
reset_timer!(to)

run_simple_evol(Run, Param, Tfinal)
println(to)
