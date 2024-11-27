push!(LOAD_PATH, "../")
using WavKinS
using GLMakie


function nk_test1(kk)
    nk = kk^2 * exp(-kk / 50.0) / 250000
    return nk
end


function nk_test(kk)
    return  nk_test1(kk)
end



## Test collision_integrals

println("---------------------------------------------------------------------")
println("Testing collisional integral")
println("")


Ms = Vector{Float64}()
errorFluxs = Vector{Float64}()
for M ∈ 2 .^ (4:9)
    kmin = 5e-3
    kmax = 5e+3

    Nk = wave_spectrum(kmin,kmax,M);
    Run = Smoluchowski(Nk; K=WavKinS.K_Smoluchowski_one, interp_scheeme=WavKinS.lin_interp);

    kk = Nk.kk
    λ = Nk.λ
  
    @. Nk.nk = nk_test.(kk);

    WavKinS.St_k!(Run)
    
    Flux = wave_spectrum(kmin,kmax,M)
    @. Flux.nk = Run.Sk.nk .* Run.ω.(kk);
    FluxNumH = integrate(Flux)

    Ene = energy(Run)
    
    println("M in k= ", M, ", Integral flux num: dH/H=", FluxNumH / Ene)
    append!(Ms,M)
    append!(errorFluxs,abs(FluxNumH / Ene))
end
fig, ax = scatter(Ms,errorFluxs,label="with_log_bins")
ylims!(1e-4, 1e0)
ax.xlabel = "M"
ax.xscale = log10
ax.ylabel = "Sum particle mass flux / Total mass"
ax.yscale = log10
ax.title = "Particle mass flux 1D error"
lines!(ax, Ms,1e2 * Ms.^-2; label=L"M^{-2}", color="red", linestyle=:dash)
axislegend(ax, position = :rt)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")

