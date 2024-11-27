push!(LOAD_PATH, "../")
using WavKinS

function nk_test1(kx)
    return exp(-abs(kx))
end

function nk_test2(kx)
    return exp(-abs(kx)^2) * kx^2
end


function nk_test(kx)
    return nk_test2(kx)
end



## Test collision_integrals

println("---------------------------------------------------------------------")
println("Testing collisional integral")
println("")

for j = 4:9
    #j=8
    M = 2^j
    kmin = 1e-2
    kmax = 1e+1

    Nk = wave_spectrum(kmin, kmax, M)
    Run = WavKinS.NLS3D(Nk)

    kk = Nk.kk
    @. Nk.nk = nk_test.(kk)

    WavKinS.St_k!(Run; compute_Sk=true, compute_γk=false, compute_ηk=false)

    Flux = wave_spectrum(kmin, kmax, M)
    @. Flux.nk = kk^4 * Run.Sk.nk
    FluxNumH = integrate_with_log_bins(Flux)

    @. Flux.nk = kk^2 * Run.Sk.nk
    FluxNumN = integrate_with_log_bins(Flux)

    AA = total_waveaction(Run)
    Ene = energy(Run)

    println("M in k= ", M, ", Integral flux num: dN/n=", FluxNumN / AA, " dH/H=", FluxNumH / Ene)
end

println("")
println("---------------------------------------------------------------------")

## 
println("---------------------------------------------------------------------")
println("Testing collisional integral with compute_γk=true, compute_ηk=true")
println("")

for j = 4:9
    #j=8
    M = 2^j
    kmin = 1e-2
    kmax = 1e+1

    Nk = wave_spectrum(kmin, kmax, M)
    Run = WavKinS.NLS3D(Nk)

    kk = Nk.kk
    @. Nk.nk = nk_test.(kk)

    WavKinS.St_k!(Run; compute_Sk=true, compute_γk=true, compute_ηk=true)

    Flux = wave_spectrum(kmin, kmax, M)
    @. Flux.nk = kk^4 * Run.Sk.nk
    FluxNumH = integrate_with_log_bins(Flux)

    @. Flux.nk = kk^2 * Run.Sk.nk
    FluxNumN = integrate_with_log_bins(Flux)

    AA = total_waveaction(Run)
    Ene = energy(Run)

    println("M in k= ", M, ", Integral flux num: dN/n=", FluxNumN / AA, " dH/H=", FluxNumH / Ene)




end

println("")
println("---------------------------------------------------------------------")
