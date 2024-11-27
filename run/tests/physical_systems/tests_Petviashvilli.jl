push!(LOAD_PATH, "../")
using WavKinS
using GLMakie
using Random


function nk_test1(kx, ky)
    return exp(-kx - abs(ky))
end

function nk_test2(kx, ky)
    return abs(ky)^3 * exp(-kx - abs(ky)) * kx^1.5 / (1.0 + abs(ky))
end


function nk_test(kx, ky)
    return nk_test2(kx, ky)
end



## Test collision_integrals

println("---------------------------------------------------------------------")
println("Testing collisional integral")
println("")


Ms = Vector{Float64}()
errorFluxEnes = Vector{Float64}()
errorFluxΩs = Vector{Float64}()
errorFluxΦs = Vector{Float64}()
for M ∈ 2 .^ (4:11)

    Mx = M
    My = M
    kxmin = 5e-3
    kxmax = 5e+3
    kymin = 5e-3
    kymax = 5e+3

    Nk = wave_spectrum_khkz(kxmin, kxmax, Mx, kymin, kymax, My)
    Run = Petviashvilli_Asymp(Nk; interp_scheeme=WavKinS.bilin_interp_khkz)


    kx = Nk.kkh
    ky = Nk.kkz
    kk = Nk.kk
    KKX = kx .* ones(length(ky))'
    KKY = ones(length(kx)) .* ky'
    @. Nk.nk = nk_test.(KKX, KKY)

    WavKinS.St_k!(Run)

    get_global_diagnostics!(Run)
    Ω = Run.diags.glob_diag["Mx"].out[end] # Total potential enstrophy
    Φ= Run.diags.glob_diag["Phi"].out[end] # Total zonostrophy
    Ene= Run.diags.glob_diag["H"].out[end] # Total energy

    integ = WavKinS.integrate_with_log_bins_khkz()
    Flux = wave_spectrum_khkz(kxmin, kxmax, Mx, kymin, kymax, My) #a field to test the convergence 

    # Check conservation of energy
    @. Flux.nk = Run.Sk.nk * WavKinS.ω_Petviashvilli_Asymp(KKX, KKY)
    Sum_Flux_Ene = integrate(integ, Flux)

     # Check conservation of pontential entrophy
    @. Flux.nk = Run.Sk.nk * WavKinS.ρ_Potential_Petviashvilli_Asymp(KKX, KKY)
    Sum_Flux_Ω = integrate(integ, Flux)

     # Check conservation of zonostrophy
    @. Flux.nk = Run.Sk.nk * WavKinS.ρ_zonostrophy_Petviashvilli_Asymp(KKX, KKY)
    Sum_Flux_Φ = integrate(integ, Flux)

    println("Mx = ", Mx, ", My = ", My, ": sum fluxes E, Ω and Φ: ", Sum_Flux_Ene / Ene, " ", Sum_Flux_Ω / Ω, " ", Sum_Flux_Φ / Φ)
    append!(Ms, M)
    append!(errorFluxEnes, abs(Sum_Flux_Ene / Ene))
    append!(errorFluxΩs, abs(Sum_Flux_Ω / Ω))
    append!(errorFluxΦs, abs(Sum_Flux_Φ / Φ))
end

## 
fig, ax = scatter(Ms, errorFluxEnes, label="error energy")
scatter!(Ms, errorFluxΩs, label="error momentum")
scatter!(Ms, errorFluxΦs, label="error zonostrophy")
ylims!(1e-5, 1e2)
ax.xlabel = "M"
ax.xscale = log10
ax.ylabel = "Error"
ax.yscale = log10
ax.title = "Conservation errors"
lines!(ax, Ms, 1e3 * Ms .^ -2; label=L"M^{-2}", color="red", linestyle=:dash)
axislegend(ax, position=:rt)
display(fig)

println("")
println("---------------------------------------------------------------------")


##
#for ix in Run.partition[1][1],iy in Run.partition[1][2]

#  println("ix =",ix," iy =",iy)
#end