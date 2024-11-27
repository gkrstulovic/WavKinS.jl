push!(LOAD_PATH, "../")
using WavKinS
using GLMakie;
GLMakie.activate!();
using TimerOutputs

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

    V0GP = 3 * sqrt(Run.c / 2.0) / 4.0

    if Run.name == "Acoustic2D"
        Cd = 4^(1 / 4) * sqrt(Run.c) / (pi * V0GP)
        Prefactor = Cd * sqrt(Run.a * Run.c)
        xP = 1
        iP = findfirst(x -> x >= 5e-2, Run.Nk.kk)
    elseif Run.name == "Acoustic3D"
        Prefactor = sqrt(3 * Run.c / (32 * V0GP^2 * π * (π + 4 * log(2) - 1.0)))
        xP = 3 / 2
        iP = findfirst(x -> x >= 5e-2, Run.Nk.kk)
    else
        Prefactor = 1.0
        xP = 0.0
    end

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

        @timeit to "Plotting" if t >= tplotEvol
            P0 = Run.diags.sp_outs["Pk"].sp[iP]
            plot_energy!(Run; fig=figE, ax=axE)
            plot_theo!(axE, kk, Prefactor * sqrt(abs(P0)) * kk .^ (-xP))
            plot_energy_flux!(Run; fig=figP, ax=axP)
            sleep(0.0001)
            tplotEvol += Param.tplot
        end


    end
end



## #########################################################################################################
##########################            Defining the simulation                     #########################
###########################################################################################################

M = 1024
kmin = 1e-3
kmax = 1e0
Nk = wave_spectrum(kmin, kmax, M)

# The Acoustic2D or Acoustic3D structures contains defines all the methods for integration, see help ?Acoustic2D ?Acoustic3D 
Run = Acoustic2D(Nk; a=1.0, interp_scheeme=WavKinS.lin_interp, time_stepping_scheeme=WavKinS.RK2_step)
kk = Nk.kk
λ = Nk.λ


kd = 0.25 * kmax
lapPower = 4
kf = 2 * kmin
Δkf = (log(kmax) - log(kmin)) / 25

@. Run.Nk.nk = exp(-0.5 * ((log(kk) - log(kf)) / Δkf)^2);
EE = energy(Run)
@. Run.FD.f = Run.Nk.nk / EE
@. Run.FD.D = (kk / kd)^(2 * lapPower)

Tfinal = 10.0; # Final time of the simulation
dt = 0.0025 # time step of the simulation


## #########################################################################################################
##########################            Outputs and plots                           #########################
###########################################################################################################

isrestart = false # set this to true if simulations uses a restart

tplot = 0.1; # Plot every tplot times
tglobal = 0.1; # Compute, store and write global quantities every tglobal times 
tspstore = 20000 # Compute and store spectral quantities every tspstore times
tspwrite = 0.1 # write spectra every tspwrite

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


## #########################################################################################################
##########################     Prepare I/0, restart, and compute diagnostics      #########################
###########################################################################################################


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


println("Run in resolultion M=", M, " with λ=", Run.Nk.λ)
println("Total wave action =", Run.diags.glob_diag["N"].out[end])
println("Total energy =", Run.diags.glob_diag["H"].out[end])


## #######################################################################################################
##########################            Make plots                                #########################
#########################################################################################################




figE, axE = plot_energy!(Run; ylims=(1e-1, 1e+3))
figP, axP = plot_energy_flux!(Run, ylims=(-0.5, 1.5))
display(figP) # figE (figP) if you want to plot energy spectrum (energy flux)


#########################################################################################################
##########################            Running the code                          #########################
#########################################################################################################


run_simple_evol(Run, Param, 0.0) # Good to run first to measure time properly
reset_timer!(to)

run_simple_evol(Run, Param, Tfinal)

#save("tutorial_energy.png", figE)
#save("tutorial_flux.png", figP)

println(to)
