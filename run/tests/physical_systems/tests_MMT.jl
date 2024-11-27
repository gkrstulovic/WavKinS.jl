push!(LOAD_PATH, "../")
using WavKinS

function nk_test1(kx)
    return exp(-abs(kx) )
end

function nk_test2(kx)
    return exp(-abs(kx)) * kx^2
end


function nk_test(kx)
    return  nk_test2(kx)
end



## Test collision_integrals

println("---------------------------------------------------------------------")
println("Testing collisional integral")
println("")

for j = 4:10
    M = 2^j
    kmin = 1e-3
    kmax = 1e+2

    Nk = wave_spectrum(kmin,kmax,M)
    Run = WavKinS.MMT(Nk; interp_scheeme=WavKinS.lin_interp);

    kk = Nk.kk
    @. Nk.nk = nk_test.(kk)

    WavKinS.St_k!(Run;compute_Sk=true, compute_γk=false, compute_ηk=false)

    
    Flux = wave_spectrum(kmin,kmax,M)
    @. Flux.nk =  Run.Sk.nk
    NFlux= integrate_with_log_bins(Flux)
    
    @. Flux.nk = Run.ω(kk) * Run.Sk.nk
    EFlux= integrate_with_log_bins(Flux)

    AA = total_waveaction(Run)
    Ene = energy(Run)

    println("M = ", M, ", Integral flux num: dN/N=", NFlux/AA, " dH/H=", EFlux/Ene)

end

println("")
println("---------------------------------------------------------------------")


println("---------------------------------------------------------------------")
println("Testing collisional integral with compute_γk=true and compute_ηk=true")
println("")

for M ∈ 2 .^ (4:10)
    kmin = 1e-3
    kmax = 1e+2

    Nk = wave_spectrum(kmin,kmax,M)
    Run = WavKinS.MMT(Nk; interp_scheeme=WavKinS.lin_interp);

    kk = Nk.kk
    @. Nk.nk = nk_test.(kk)

    WavKinS.St_k!(Run;compute_Sk=true, compute_γk=true, compute_ηk=true)

    
    Flux = wave_spectrum(kmin,kmax,M)
    @. Flux.nk =  Run.Sk.nk
    NFlux= integrate_with_log_bins(Flux)
    
    @. Flux.nk = Run.ω(kk) * Run.Sk.nk
    EFlux= integrate_with_log_bins(Flux)

    AA = total_waveaction(Run)
    Ene = energy(Run)

    println("M = ", M, ", Integral flux num: dN/N=", NFlux/AA, " dH/H=", EFlux/Ene)

end

println("")
println("---------------------------------------------------------------------")