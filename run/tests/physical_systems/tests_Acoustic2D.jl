# Tests for integration and interpolation for isotropic systems

push!(LOAD_PATH, "../")
using WavKinS
using GLMakie
using Random

function nk_test(kk)
    nk = kk^2 * exp(-kk / 50.0) / 250000
    return nk
end

function Stk_test(kk)# Stk_test=Stk1_test +Stk2_test
    Sfactor = 2 * pi / (sqrt(6) * 0.5 * 1)
    V0GP = 3 * sqrt(1 / 2.0) / 4.0
    Stk = 1125.0 / 8.0 + 45 * kk / 16.0 + 9 * kk^2 / 400.0 - 93.0 * kk^3 / 40000.0 + kk^7 / 8750000000000.0
    Stk = V0GP^2 * Sfactor * Stk * exp(-kk / 50.0)
    return Stk
end

function Stk1_test(k)  # Stk1_test: term with integral from 0 to k
    Sfactor = 2 * pi / (sqrt(6) * 0.5 * 1)
    V0GP = 3 * sqrt(1 / 2.0) / 4.0
    Stk1 = -(6 / 25) * exp(-k / 25) * k^2 + 6 / 25.0 * exp(-k / 50) * k^2 - 9 * exp(-k / 25) * k^3 / 2500.0 - 3 * exp(-k / 50) * k^3 / 2500.0 - 3 * exp(-k / 25) * k^4 / 125000.0 - exp(-k / 25) * k^5 / 12500000.0 + exp(-k / 50) * k^7 / 8750000000000.0
    return V0GP^2 * Sfactor * Stk1
end

function Stk2_test(k) # Stk1_test: term with integral from k to \infty
    Sfactor = 2 * pi / (sqrt(6) * 0.5 * 1)
    V0GP = 3 * sqrt(1 / 2.0) / 4.0
    Stk2 = (1125 * exp(-k / 50)) / 8 + 45 / 16 * exp(-k / 50) * k + 6 / 25 * exp(-k / 25) * k^2 - 87 / 400.0 * exp(-k / 50) * k^2 + (9 * exp(-k / 25) * k^3) / 2500.0 - (9 * exp(-k / 50) * k^3) / 8000.0 + (3 * exp(-k / 25) * k^4) / 125000.0 + (exp(-k / 25) * k^5) / 12500000.0
    return V0GP^2 * Sfactor * Stk2
end

function dStk2_test(k1, k)
    Sfactor = 2 * pi / (sqrt(6) * 0.5 * 1)
    V0GP = 3 * sqrt(1 / 2.0) / 4.0
    k2 = k1 - k
    dS = k1 * k2 * (-2 * nk_test(k2) * nk_test(k) + 2 * nk_test(k1) * nk_test(k) + 2 * nk_test(k1) * nk_test(k2))
    return dS * V0GP^2 * Sfactor
end

function dEnergy(q)
    dE = (375 * q^3) / 8.0 + (15 * q^4) / 16.0 + (33 * q^5) / 4000.0 - (9 * q^6) / 25000.0 - (9 * q^7) / 8750000.0 - (9 * q^8) / 3500000000.0 - q^9 / 175000000000.0
    return dE * exp(-q / 50.0)
end


## Test collision_integrals

println("---------------------------------------------------------------------")
println("Error computing collisional integral")
println("")

Ms = Vector{Float64}()
errorFluxs = Vector{Float64}()
for j = 4:10
    M = 2^j
    kmin = 1e-2
    kmax = 1e+4

    Nk = wave_spectrum(kmin,kmax,M)
    kk = Nk.kk
    λ = Nk.λ
    @. Nk.nk = nk_test.(kk)
    Run = Acoustic2D(Nk; interp_scheeme=WavKinS.powGauss_interp)

    WavKinS.St_k!(Run)

    StkTheo = similar(Run.Sk.nk)
    StkTheo = Stk_test.(kk)

    ErrMaxSt = maximum(abs.(StkTheo - Run.Sk.nk))

    Flux = wave_spectrum(kmin,kmax,M)
    @. Flux.nk = Flux.kk^2 * Run.Sk.nk
    FluxNum = integrate_with_log_bins(Flux)
    @. Flux.nk = Flux.kk^2 * StkTheo
    FluxNumTheo = integrate_with_log_bins(Flux)

    println("M = ", M, ": max error Stk: ", ErrMaxSt, ", Integral flux num=", FluxNum, ", Integral flux theo=", FluxNumTheo)

    fig, ax = plot_1D_base!(Run.Sk.kk, abs.((Run.Sk.nk - StkTheo)./Run.Sk.nk), xlabel=L"$k$", ylabel="Relative error", title="Collision integral error", yscale=log10)
    scatter!(ax, Run.Sk.kk, abs.((Run.Sk.nk - StkTheo)./Run.Sk.nk))
    display(fig)
    sleep(5)
    append!(Ms,M)
    append!(errorFluxs,abs(FluxNum))
end
fig, ax = scatter(Ms,errorFluxs,label="with_log_bins")
ylims!(1e3, 1e10)
ax.xlabel = "M"
ax.xscale = log10
ax.ylabel = "Sum energy flux"
ax.yscale = log10
ax.title = "Energy flux 1D error"
lines!(ax, Ms,1e10 * Ms.^-2; label=L"M^{-2}", color="red", linestyle=:dash)
axislegend(ax, position = :rt)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")



#benchmarking
#using BenchmarkTools
#@btime
