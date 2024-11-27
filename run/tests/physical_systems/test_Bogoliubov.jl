push!(LOAD_PATH, "../")
using WavKinS

function nk_test1(kx)
    return exp(-abs(kx) )
end

function nk_test2(kx)
    return exp(-abs(kx)^1) * kx^4 /(1)
end

function nk_test3(kx)
    return 1/(1 + kx^3)
end



function nk_test(kx)
    return  nk_test3(kx)
end



## Test collision_integrals

println("---------------------------------------------------------------------")
println("Testing collisional integral")
println("")

for j = 4:12
#j=12
    M = 2^j
    kmin = 1e-1
    kmax = 1e+2
    ξn = (1/1.e-0)

    Nk = wave_spectrum(kmin,kmax,M)
    Run = WavKinS.Bogoliubov3D(Nk; ξ = ξn,interp_scheeme=WavKinS.lin_interp);
    #Run = Acoustic3D(Nk; interp_scheeme=WavKinS.powGauss_interp)


    kk = Nk.kk
    @. Nk.nk = nk_test.(kk)

    WavKinS.St_k!(Run)
    
    Flux = wave_spectrum(kmin,kmax,M)
    @. Flux.nk = WavKinS.ω_Bogoliubov(kk; ξ = ξn) .* kk.^2 * Run.Sk.nk
    FluxNum= integrate_with_log_bins(Flux)

    println("M = ", M, ", Integral flux num=", FluxNum)

end

println("")
println("---------------------------------------------------------------------")


##
#for ix in Run.partition[1][1],iy in Run.partition[1][2]

  #  println("ix =",ix," iy =",iy)
#end