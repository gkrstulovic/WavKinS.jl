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

    @time while t < tini + Tfinal
      #  println(" Min nk=",minimum(Run.Nk.nk)," Max nk=",maximum(Run.Nk.nk))
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


###########################################################################################################
###########################################################################################################
## Defining the simulation

Mx = 128
My = 128
kxmin = 1e-4
kxmax = 1.e-0
kymin = 1e-4
kymax = 1e-0

Nk = wave_spectrum_khkz(kxmin,kxmax,Mx,kymin,kymax,My);

# The Petviashvilli_Asymp structure contains all the methods for integration, see help ?Petviashvilli_Asymp for options
Run = Petviashvilli_Asymp(Nk; interp_scheeme=WavKinS.bilin_interp_khkz, drive_scheeme=WavKinS.AB_Euler_step);


###########################################################################################################
##########################            running parmeters                           #########################
###########################################################################################################

isrestart = false # set this to true if simulations uses a restart

Tfinal = 10.0; # Final time of the simulation
dt = 0.025; # time step of the simulation

tplot = 0.1; # Plot every tplot times
tglobal = 0.00025; # Compute, store and write global quantities every tglobal times 
tspstore = 20000 # Compute and store spectral quantities every tspstore times

tspwrite = 0.1# write spectra every tspwrite
outputDir = "./" #output directory
write_global = true #write global
write_nk = true #write wave action
write_Ek = true #write energy spectrum
write_Pk = false #write energy flux   

# Creating simulations parameters and setting write_sp for each spectrum
write_spectral = write_nk || write_Ek || write_Pk
Param = simulation_parameters(dt, tplot, tglobal, tspstore, tspwrite,
    outputDir, write_global, write_spectral)

Run.diags.sp_outs["nk"].write_sp = write_nk
Run.diags.sp_outs["Ek"].write_sp = write_Ek
Run.diags.sp_outs["Pkh"].write_sp = write_Pk
Run.diags.sp_outs["Pkz"].write_sp = write_Pk

###########################################################################################################
###############     Setting forcing, dissipation and initial condition            #########################
###########################################################################################################


kkx = Nk.kkh
kky = Nk.kkz
KX = kkx .* ones(length(kky))';
KY = ones(length(kkx)) .* kky';
K2 = KX .^ 2 .+ KY .^ 2;

nu = 1
nuinv = 1
kxinf = 1*kxmin
kxsup = .2*kxmax 
kyinf =  5*kymin
kysup =  .2*kymax 
lapPower = 2
lapInvPower = 3


kxf =5e-2;
kyf =5e-2;
Δkfy = 0.005;
Δkfx = 0.005;


@. Run.Nk.nk = exp.(-.5 * ((KX - kxf)^2/Δkfx^2 +  (KY - kyf)^2/Δkfy^2) );# * exp(- (KY/kysup)^2);
AA = total_waveaction(Run)
@. Run.Nk.nk = Run.Nk.nk / AA;
@. Run.FD.f = Run.Nk.nk;
#@. Run.Nk.nk += 1e-80;
@. Run.FD.D = nu * (((KX / kxsup)^2 + (KY / kysup)^2)^(lapPower));
@. Run.FD.D += nuinv * (  (KX / kxinf) .^ (-2 * lapInvPower) + (KY / kyinf) .^ (-2 * lapInvPower));



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
println("Run in resolultion Mx=", Mx, " with λx=", Nk.logλh)
println("                   My=", My, " with λx=", Nk.logλz)
println("Total wave action =", Run.diags.glob_diag["N"].out[end])
println("Total energy =", Run.diags.glob_diag["H"].out[end])


fig, ax, hm = plot_energy!(Run; zlims=(1e-8, 1e3))
display(fig)


## ## ####################################################################################################
##########################            Running the code                          #########################
#########################################################################################################

run_simple_evol(Run, Param, 0.0) # Good to run first to measure time properly
reset_timer!(to)

run_simple_evol(Run, Param, Tfinal)
println(to)

H = Run.diags.glob_diag["H"].out;
ΔH = (H[end]-H[1])/H[1]
println("Energy conservation =",ΔH)
