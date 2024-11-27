push!(LOAD_PATH, "../../")
using WavKinS
using TimerOutputs
using GLMakie;
GLMakie.activate!();

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

    dtmin = 1.e-3
    dtmax = 0.95 / maximum(0.25 * Run.FD.D)

    V0GP = 3 * sqrt(Run.c / 2.0) / 4.0
    d = Run.dimension
    dΩ = Run.dΩ
    iP = findfirst(x -> x >= 1.e-0, Run.Nk.kk)



    PrefactorOLD = sqrt(3 * Run.c / (32 * V0GP^2 * π * (π + 4 * log(2) - 1.0)))
    dIkzZS = 64 * (-1 + π + log(16)) / 3
    dIkaco = 105.149
    Prefactor = sqrt(2 / (π * V0GP^2 * Run.c * dIkaco))
    println(PrefactorOLD, " ", Prefactor)
    xP = 3 / 2

    dIkz = 0.738008
    Prefactor2 = sqrt(9 * Run.c^2 * Run.ξ^2 / (4 * sqrt(2) * π * dIkz * V0GP^2))

    @time while t < tini + Tfinal
        Param.dt = adaptative_time_step(Run, dtmin, dtmax, Param.dt)
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

        @timeit to "Storing spectral quantities" if t >= tsptoreEvol
            store_spectral!(Run)
            tsptoreEvol += Param.tspstore
        end

        @timeit to "Plotting" if t >= tplotEvol
            if size(Run.diags.sp_store["Pk"], 2) >= 1
                P0 = Run.diags.sp_store["Pk"][iP, end]
            else
                P0 = 1e-2
            end
            #P0 = Run.diags.glob_diag["Disp"].out[end]
            plot_energy!(Run; fig=figE, ax=axE)
            plot_theo!(axE, kk, Prefactor * sqrt(abs(P0)) * kk .^ (-xP))
            plot_theo!(axE, kk, Prefactor2 * sqrt(abs(P0)) * kk .^ (1), color="green")
            plot_energy_flux!(Run; fig=figP, ax=axP)
            sleep(0.0001)
            tplotEvol += Param.tplot
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


M = 128
kmin = 1e-3
kmax = 1.

Nk = wave_spectrum(kmin, kmax, M)

# The Bogoliubov3D structure contains all the methods for integration, see help ?Bogoliubov3D for options
Run = Bogoliubov3D(Nk; ξ=(1 / 1), interp_scheeme=WavKinS.lin_interp, time_stepping_scheeme=WavKinS.RK2_step)


###########################################################################################################
##########################            running parmeters                           #########################
###########################################################################################################

isrestart = false # set this to true if simulations uses a restart

Tfinal = 20; # Final time of the simulation
dt = 0.01 # time step of the simulation2
tplot = 0.5; # Plot every tplot times
tglobal = 0.5; # Compute, store and write global quantities every tglobal times 
tspstore = 0.5 # Compute and store spectral quantities every tspstore times

tspwrite = 1 # write spectra every tspwrite
outputDir = "./" #output directory
write_global = true #write global
write_nk = true #write wave action
write_Ek = true #write energy spectrum
write_Pk = true #write energy flux   

# Creating simulations parameters and setting write_sp for each spectrum
write_spectral = write_nk || write_Ek || write_Pk
Param = simulation_parameters(dt, tplot, tglobal, tspstore, tspwrite,
    outputDir, write_global, write_spectral)

Run.diags.sp_outs["nk"].write_sp = write_nk
Run.diags.sp_outs["Ek"].write_sp = write_Ek
Run.diags.sp_outs["Pk"].write_sp = write_Pk

###########################################################################################################
###############     Setting forcing, dissipation and initial condition            #########################
###########################################################################################################
#

kk = Nk.kk
nu = 1
kd = 0.1 * kmax
lapPower = 2
LSInvlapPower = 1

kf = 1.e-4
Δkf = 2e-3

@. Run.Nk.nk = kk^2 * exp(-0.5 * (kk - kf)^2 / (Δkf)^2);

AA = total_waveaction(Run)
@. Run.Nk.nk = Run.Nk.nk / AA
EE = energy(Run)
@. Run.FD.f = Run.Nk.nk / EE
@. Run.FD.D = nu * (kk / kd)^(2 * lapPower)

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


figE, axE = plot_energy!(Run; ylims=(1e-2, 1e4))
figP, axP = plot_energy_flux!(Run, ylims=(-1.5, 1.5))
display(figE) # figE (figP) if you want to plot energy spectrum (energy flux)


## ######################################################################################################
##########################            Running the code                          #########################
#########################################################################################################

run_simple_evol(Run, Param, 0.0) # Good to run first to measure time properly
reset_timer!(to)

run_simple_evol(Run, Param, Tfinal)
println(to)

