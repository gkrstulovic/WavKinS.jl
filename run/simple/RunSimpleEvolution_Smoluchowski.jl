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

    dtmin = 1.e-3
    dtmax = 1 / maximum( .25*Run.FD.D)

    d = Run.dimension
    dΩ = Run.dΩ

    #if Run.ω == K_Smoluchowski_multiplicative
    Prefactor = 1/sqrt(2*pi)
    xP = 3/2
    iP = findfirst(x -> x >= 5e-2, Run.Nk.kk)


    @time while t < tini + Tfinal
        if any(isnan.(Run.Nk.nk))
            @warn "We got NaNs. Stop the temporal loop"
            break
        end


        @timeit to "Advance" advance!(step_scheeme, Run, Param.dt)
        t = Run.t
        #   Param.dt = adaptative_time_step(Run, dtmin, dtmax, Param.dt)

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
kmin = 5e-3
kmax = 1e0
Nk = wave_spectrum(kmin, kmax, M)

# The Smoluchowski structure contains defines all the methods for integration, see help ?Smoluchowski
Run = Smoluchowski(Nk; K=WavKinS.K_Smoluchowski_multiplicative,interp_scheeme=WavKinS.lin_interp, time_stepping_scheeme=WavKinS.RK2_step)
kk = Nk.kk
λ = Nk.λ


kd = 0.25 * kmax
lapPower = 4
kf = 3 * kmin
Δkf = (log(kmax) - log(kmin)) / 25

@. Run.Nk.nk = kk^2 * exp(-0.5 * ((log(kk) - log(2 * kmin)) / Δkf)^2);
EE = energy(Run)
@. Run.FD.f = Run.Nk.nk / EE
@. Run.FD.D = (kk / kd)^(2 * lapPower)

Tfinal = 50.0; # Final time of the simulation
dt = 0.025 # time step of the simulation


## #########################################################################################################
##########################            Outputs and plots                           #########################
###########################################################################################################

isrestart = false # set this to true if simulations uses a restart

tplot = 1.0; # Plot every tplot times
tglobal = 1.0; # Compute, store and write global quantities every tglobal times 
tspstore = 20000 # Compute and store spectral quantities every tspstore times
tspwrite = 1.0 # write spectra every tspwrite

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

figE, axE = plot_energy!(Run; ylims=(1e-2, 1e+4))
figP, axP = plot_energy_flux!(Run, ylims=(-0.5, 1.5))
display(figE) # figE (figP) if you want to plot energy spectrum (energy flux)


#########################################################################################################
##########################            Running the code                          #########################
#########################################################################################################

run_simple_evol(Run, Param, 0.0) # Good to run first to measure time properly
reset_timer!(to)

run_simple_evol(Run, Param, Tfinal)
println(to)

#save("Smoluchowski_spectrum.png", figE)