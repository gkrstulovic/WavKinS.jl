push!(LOAD_PATH, "../")
using WavKinS
using GLMakie; GLMakie.activate!()
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

    dtmin = 1.e-6
    dtmax = 1 / maximum( .25*Run.FD.D)

    @time while t < tini + Tfinal
    #  println(" Min nk=",minimum(Run.Nk.nk)," Max nk=",maximum(Run.Nk.nk))
        @timeit to "Advance" advance!(step_scheeme, Run, Param.dt)
        t = Run.t
        Param.dt = adaptative_time_step(Run, dtmin, dtmax, Param.dt)

        @timeit to "Global quantities" if t >= tglobalEvol
            get_global_diagnostics!(Run)
            if Param.write_global
                output_global(Run, Param.outputDir)
            end
            Energy = Run.diags.glob_diag["H"].out[end]
            WaveAction = Run.diags.glob_diag["N"].out[end]
            conservation_ratio = energy_conservation_ratio(Run)     
            println("t =", Run.t, "   H(t) =", Energy, "     N(t) = ", WaveAction, "     Energy conservation ratio = ", conservation_ratio)
            tglobalEvol += Param.tglobal
        end

        @timeit to "Storing spectral quantities" if t >= tsptoreEvol
            store_spectral!(Run)
            tsptoreEvol += Param.tspstore
        end

        @timeit to "Plotting" if t >= tplotEvol
            plot_energy!(Run; fig=fig, ax=ax, hm=hm)
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


## #########################################################################################################
##########################            Defining the simulation                     #########################
###########################################################################################################

Mh = 16
Mz = 16
khmin = 5e-3
khmax = 1e0
kzmin = 5e-3
kzmax = 1e0
Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz);

# The Stratified_Asymp structure contains all the methods for integration, see help ?Stratified_Asymp for options
Run = Stratified_Asymp(Nk; interp_scheeme=WavKinS.bilin_interp_khkz, time_stepping_scheeme=WavKinS.RK2_step);
kkh = Nk.kkh;
kkz = Nk.kkz;
kk = Nk.kk;
λh = Nk.λh;
λz = Nk.λz;
KH = kkh .* ones(length(kkz))';
KZ = ones(length(kkh)) .* kkz';
K2 = KH .^ 2 .+ KZ .^ 2;


kdsup = 0.5*khmax
kdinf = 1*khmin 
lapPower = 4
lapInvPower = lapPower
of = 1.0
kf = 0.03
thf = atan(of);
khf = kf*sin(thf);
kzf = kf*cos(thf);
Δkf = (log(khmax) - log(khmin)) / 50;

@. Run.Nk.nk = exp.(-((log.(KH) - log(khf))^2 + (log.(KZ) - log(kzf))^2) / Δkf^2) ./ Run.ω.(KH,KZ) # localized
EE = energy(Run)
@. Run.FD.f = Run.Nk.nk / EE 
@. Run.FD.D = ( (K2/kdsup^2).^(lapPower) +  (KH / kdinf) .^ (-2 * lapInvPower) + (KZ / kdinf) .^ (-2 * lapInvPower) ) ./ Run.ω.(KH,KZ)

Tfinal = 1.0; # Final time of the simulation
dt = 1e-6; # time step of the simulation


## #########################################################################################################
##########################            Outputs and plots                           #########################
###########################################################################################################

isrestart = false # set this to true if simulations uses a restart

tplot = 0.1; # Plot every tplot times
tglobal = 0.01; # Compute, store and write global quantities every tglobal times 
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
Run.diags.sp_outs["Pkh"].write_sp = write_Pk
Run.diags.sp_outs["Pkz"].write_sp = write_Pk 


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


println("Run in resolultion Mh=", Mh, " with λh=", λh)
println("                   Mz=", Mz, " with λz=", λz)
println("Total wave action =", Run.diags.glob_diag["N"].out[end])
println("Total energy =", Run.diags.glob_diag["H"].out[end])


## #######################################################################################################
##########################            Make plots                                #########################
#########################################################################################################

fig, ax, hm = plot_energy!(Run; zlims=(1e-6,1e3))
display(fig)


#########################################################################################################
##########################            Running the code                          #########################
#########################################################################################################

run_simple_evol(Run, Param, 0.0) # Good to run first to measure time properly
reset_timer!(to)

run_simple_evol(Run, Param, Tfinal)
println(to)



