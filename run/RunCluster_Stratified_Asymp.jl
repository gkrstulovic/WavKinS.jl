push!(LOAD_PATH, "../")
using WavKinS
using TimerOutputs

## Init basics
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

        @timeit to "Advance" advance!(step_scheeme, Run, Param.dt)
        t = Run.t
        Param.dt = adaptative_time_step(Run, dtmin, dtmax, Param.dt)

        @timeit to "Global quantities" if t >= tglobalEvol
            AA = total_waveaction(Run)
            Energy = energy(Run)
            Disp = energy_dissipation(Run)
            
            push!(Run.diags.glob_diag["Times"].out, t)
            push!(Run.diags.glob_diag["N"].out, AA)
            push!(Run.diags.glob_diag["H"].out, Energy)
            push!(Run.diags.glob_diag["Disp"].out, Disp)
            
            if Param.write_global
                output_global(Run, Param.outputDir)
            end
            tglobalEvol += Param.tglobal

            # Check energy conservation
            conservation_ratio = energy_conservation_ratio(Run)     
            println("t =", t, "   H(t) =", Energy, "     N(t) = ", AA, "     Energy conservation ratio = ", conservation_ratio)
        end

        @timeit to "Storing spectral quantities" if t >= tsptoreEvol
            Run.diags.sp_store["nk"] = cat(Run.diags.sp_store["nk"], Run.Nk.nk; dims=3)
            Run.diags.sp_store["Pkh"] = cat(Run.diags.sp_store["Pkh"], zeros(Run.Nk.Mh,Run.Nk.Mz); dims=3)
            Run.diags.sp_store["Pkz"] = cat(Run.diags.sp_store["Pkz"], zeros(Run.Nk.Mh,Run.Nk.Mz); dims=3)
            Pkh = Run.diags.sp_store["Pkh"][:,:,end]
            Pkz = Run.diags.sp_store["Pkz"][:,:,end]
            energy_flux!([Pkh, Pkz], Run)
          
            tsptoreEvol += Param.tspstore
        end

        @timeit to "Writing spectral quantities" if t >= tspwriteEvol
            if Param.write_spectral
                if Run.diags.sp_outs["nk"].write_sp
                    Run.diags.sp_outs["nk"].sp .= Run.Nk.nk
                end
                if Run.diags.sp_outs["Ek"].write_sp
                    energy_spectrum!(Run.diags.sp_outs["Ek"].sp, Run)
                end
                if Run.diags.sp_outs["Pkh"].write_sp
                    Pkh = Run.diags.sp_outs["Pkh"].sp
                    Pkz = Run.diags.sp_outs["Pkz"].sp
                    energy_flux!([Pkh, Pkz], Run)
                end
                output_spectra(Run, Param.outputDir)
            end
            tspwriteEvol += Param.tspwrite
        end
    end
end


###########################################################################################################
## Defining the simulation
###########################################################################################################

# Parsing arguments
sim_type = ARGS[1]
println(sim_type)
M = parse(Int64, ARGS[2])

if sim_type == "big"
    kmin = 1e-4
    kmax = 1e0
else
    kmin = 5e-3
    kmax = 1e0
end
   
kdinf = parse(Float64, ARGS[3])
kdsup = parse(Float64, ARGS[4]) 


###########################################################################################################
##########################            running parmeters                           #########################
###########################################################################################################

isrestart = true # set this to true if simulations uses a restart

Mh = M
Mz = M
khmin = kmin
khmax = kmax
kzmin = kmin
kzmax = kmax
lapInvPower = 4
lapPower = lapInvPower   

Tfinal = 500; # Final time of the simulation
dt = 1e-6; # time step of the simulation

tglobal = 1e-2 # Compute, store and write global quantities every tglobal times 
tspstore = 2e-2 # Compute and store spectral quantities every tspstore times
tspwrite = 2e-2 # write spectra every tspwrite

# output directory and other parameters
if sim_type == "localized" || sim_type == "big" # Anisotropic localized forcing (also used for big simulation)
    kfh = parse(Float64, ARGS[5]) 
    kfz = parse(Float64, ARGS[6]) 
    outputDir = "/scratch/vlabarre/WKE_Stratified_Asymp/$(sim_type)_M$(M)_kmin$(khmin)_kmax$(khmax)_kdinf$(kdinf)_kdsup$(kdsup)_kfh$(kfh)_kfz$(kfz)_lapPower$(lapPower)/" 
elseif  sim_type == "isotropic" # Isotropic forcing
    kf = parse(Float64, ARGS[5]) 
    outputDir = "/scratch/vlabarre/WKE_Stratified_Asymp/$(sim_type)_M$(M)_kmin$(khmin)_kmax$(khmax)_kdinf$(kdinf)_kdsup$(kdsup)_kf$(kf)_lapPower$(lapPower)/"
elseif sim_type == "decaying"   # Initial log normal
    outputDir = "/scratch/vlabarre/WKE_Stratified_Asymp/$(sim_type)_M$(M)_kmin$(khmin)_kmax$(khmax)_kdinf$(kdinf)_kdsup$(kdsup)_lapPower$(lapPower)/"
elseif sim_type == "decayingKZ"   # Initial KZ spectrum
    outputDir = "/scratch/vlabarre/WKE_Stratified_Asymp/$(sim_type)_M$(M)_kmin$(khmin)_kmax$(khmax)_kdinf$(kdinf)_kdsup$(kdsup)_lapPower$(lapPower)/"
elseif sim_type == "decayingLD"   # Initial Lvov-Dematteis spectrum
    outputDir = "/scratch/vlabarre/WKE_Stratified_Asymp/$(sim_type)_M$(M)_kmin$(khmin)_kmax$(khmax)_kdinf$(kdinf)_kdsup$(kdsup)_lapPower$(lapPower)/"
elseif sim_type == "decayingLC"   # Initial Lanchon-Cortet spectrum
    outputDir = "/scratch/vlabarre/WKE_Stratified_Asymp/$(sim_type)_M$(M)_kmin$(khmin)_kmax$(khmax)_kdinf$(kdinf)_kdsup$(kdsup)_lapPower$(lapPower)/"
elseif sim_type == "decayingGM"   # Initial Garrett-Munk spectrum
    outputDir = "/scratch/vlabarre/WKE_Stratified_Asymp/$(sim_type)_M$(M)_kmin$(khmin)_kmax$(khmax)_kdinf$(kdinf)_kdsup$(kdsup)_lapPower$(lapPower)/"
else
    @error "Unvalid sim_type (first argument of the script)" 
end


write_global = true #write global
write_nk = true #write wave action
write_Ek = true #write energy spectrum
write_Pk = false #write energy flux   
write_spectral = write_nk || write_Ek || write_Pk

Param = simulation_parameters(dt, 0.0, tglobal, tspstore, tspwrite,
    outputDir, write_global, write_spectral)

Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz);
Run = Stratified_Asymp(Nk; interp_scheeme=WavKinS.bilin_interp_khkz, time_stepping_scheeme=WavKinS.RK2_step);

Run.diags.sp_outs["nk"].write_sp = write_nk
Run.diags.sp_outs["Ek"].write_sp = write_Ek
Run.diags.sp_outs["Pkh"].write_sp = write_Pk
Run.diags.sp_outs["Pkz"].write_sp = write_Pk 
    

###########################################################################################################
###############     Setting forcing, dissipation and initial condition            #########################
###########################################################################################################

kkh = Nk.kkh;
kkz = Nk.kkz;
kk = Nk.kk;
λh = Nk.λh;
λz = Nk.λz;
KH = kkh .* ones(length(kkz))';
KZ = ones(length(kkh)) .* kkz';
K2 = KH .^ 2 .+ KZ .^ 2;


if sim_type == "localized" # Anisotropic localized forcing
    #Δkf = 1.0*khmin;
    #@. Run.Nk.nk = exp.(-((KH - kfh)^2 + (KZ - kfz)^2) / (Δkf)^2) ./ Run.ω.(KH,KZ)
    Δkf = exp((log(khmax) - log(khmin))/50);
    @. Run.Nk.nk = exp.(-((log.(KH) - log(kfh))^2 + (log.(KZ) - log(kfz))^2) / (log(Δkf))^2) ./ Run.ω.(KH,KZ) # localized
    EE = energy(Run);
    @. Run.FD.f = Run.Nk.nk / EE;
    @. Run.Nk.nk = 0.0;

elseif  sim_type == "isotropic" # Isotropic forcing
    Δkf = 1.0*khmin;
    @. Run.Nk.nk = exp(-(sqrt(K2) - kf)^2 / (Δkf)^2) ./ Run.ω.(KH,KZ)
    EE = energy(Run);
    @. Run.FD.f = Run.Nk.nk / EE;
    @. Run.Nk.nk = 0.0 ;

elseif sim_type == "decaying"   # Initial log normal
    @. Run.Nk.nk = exp.(-(log10.(K2)/2+1.0).^2/0.005) ./ Run.ω.(KH,KZ)
    EE = energy(Run);
    @. Run.Nk.nk = Run.Nk.nk / EE;
    @. Run.FD.f = 0.0;

elseif sim_type == "decayingKZ"   # Initial KZ spectrum
    @. Run.Nk.nk = KH.^(-3.5) .* KZ.^(-0.5)
    EE = energy(Run);
    @. Run.Nk.nk = Run.Nk.nk / EE;
    @. Run.FD.f = 0.0;

elseif sim_type == "decayingLD"   # Initial Lvov-Dematteis spectrum
    @. Run.Nk.nk = KH.^(-3.69)
    EE = energy(Run);
    @. Run.Nk.nk = Run.Nk.nk / EE;
    @. Run.FD.f = 0.0;

elseif sim_type == "decayingLC"   # Initial Lanchon-Cortet spectrum
    @. Run.Nk.nk = KH.^(-3) .* KZ.^(-1)
    EE = energy(Run);
    @. Run.Nk.nk = Run.Nk.nk / EE;
    @. Run.FD.f = 0.0;

elseif sim_type == "decayingGM"   # Initial Garrett-Munk spectrum
    @. Run.Nk.nk = KH.^(-4) .* KZ.^(0)
    EE = energy(Run);
    @. Run.Nk.nk = Run.Nk.nk / EE;
    @. Run.FD.f = 0.0;

elseif sim_type == "big" # Anisotropic localized forcing
    Δkf = exp((log(khmax) - log(5e-3))/50);
    @. Run.Nk.nk = exp.(-((log.(KH) - log(kfh))^2 + (log.(KZ) - log(kfz))^2) / (log(Δkf))^2) ./ Run.ω.(KH,KZ) # localized
    EE = energy(Run);
    @. Run.FD.f = Run.Nk.nk / EE;
    @. Run.Nk.nk = 0.0;

else
    @error "Unvalid sim_type (first argument of the script)" 
end


# Dissipation 
@. Run.FD.D = ( (K2/kdsup^2).^(lapPower) +  (KH / kdinf) .^ (-2 * lapInvPower) + (KZ / kdinf) .^ (-2 * lapInvPower) ) ./ Run.ω.(KH,KZ);

Nini = total_waveaction(Run);
Hini = energy(Run);
Disp = energy_dissipation(Run);


println("Run in resolultion Mh=", Mh, " with λh=", λh)
println("                   Mz=", Mz, " with λz=", λz)
println("Total wave action =", Nini)
println("Total energy =", Hini)

if ~isrestart || ~isdir(Param.outputDir)
    mkdir(Param.outputDir)

    push!(Run.diags.glob_diag["Times"].out, 0.0)
    push!(Run.diags.glob_diag["N"].out, Nini)
    push!(Run.diags.glob_diag["H"].out, Hini)
    push!(Run.diags.glob_diag["Disp"].out, Disp)
    
    init_IO(Run, Param.outputDir)
else
    load_spectrum(Run, Param.outputDir)
end


## ## ####################################################################################################
##########################            Running the code                          #########################
#########################################################################################################

run_simple_evol(Run, Param, 0.0) # Good to run first to measure time properly
reset_timer!(to)

run_simple_evol(Run, Param, Tfinal)
println(to)
