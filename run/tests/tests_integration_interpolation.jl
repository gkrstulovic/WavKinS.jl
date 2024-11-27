# Tests for integration and interpolation for isotropic systems

push!(LOAD_PATH, "../")
using WavKinS
using GLMakie
using Random

function nk_test(kk)
    nk = kk^2 * exp(-kk / 50.0) / 250000
    return nk
end

function nk_test2(kk)
    nk = kk^2 * exp(-(kk - 0.1)^2 / 0.01^2)
    return nk
end

function nk_test3(kk)
    nk = kk^2 / (0.01 + kk^3)
    return nk
end

function nk_test4(kk)
    nk = 3*kk^2/2 + kk^3 - cos(kk*pi)
    return nk
end

##   Test interpolation
println("---------------------------------------------------------------------")
println("Relative error interpolation")
println("") 

function ntest(k)
    return nk_test(k)
end

Ms = Vector{Float64}()
errorInterpIns = Vector{Float64}()
errorInterpOuts = Vector{Float64}()
for j = 4:16
    M = 2^j
    kmin = 1e-2
    kmax = 1e+3

    Nk = wave_spectrum(kmin,kmax,M)
    kk = Nk.kk
    λ = Nk.λ
    @. Nk.nk = ntest.(kk)

    interp = WavKinS.lin_interp(kk)
    #interp = WavKinS.linlog_interp(kk)
    #interp = WavKinS.powexp_interp(kk) 
    #interp = WavKinS.powGauss_interp(kk)
    #interp = WavKinS.BS_interp(kk)
    WavKinS.update_coeff_interp!(interp, Nk)

    #ktest = [2.e-3, 1.e-2, 0.12, 1.1, 4.56e+1, 23.0, 151.0, 900.0, 5000.0]
    InterpolationPoints = 500
    Random.seed!(1234)
    Irand = sort(rand(InterpolationPoints) .* (M-1) .+ 1)
    ktestIn = kk[1] * λ .^ (Irand .- 1)
    ktestOut = [kmin/2]
    NtestIn = similar(ktestIn)
    NtestOut = similar(ktestOut)

    for i in eachindex(ktestIn)
        NtestIn[i] = val_nk(interp, Nk, ktestIn[i])
    end
    for i in eachindex(ktestOut)
        NtestOut[i] = val_nk(interp, Nk, ktestOut[i])
    end

    inds = NtestIn .> 1.e-40
    errorInterpIn = abs.((NtestIn[inds] .- ntest.(ktestIn[inds])) ./ ntest.(ktestIn[inds]))
    errorInterpOut = abs.((NtestOut .- ntest.(ktestOut)) ./ ntest.(ktestOut))

    println("M = ", M, ": relative error interpolation: ", maximum(errorInterpIn), ", relative error extrapolation: ", maximum(errorInterpOut))

    fig, ax = plot_1D_base!(Nk.kk, Nk.nk; xlabel=L"$k$", ylabel=L"$n_k$", title="Interpolation 1D", xlims=(kmin / 10, 10 * kmax), ylims=(1e-20, 1e+0))
    scatter!(ax, kk, Nk.nk)
    scatter!(ax, ktestIn, NtestIn, label="Interpolation")
    scatter!(ax, ktestOut, NtestOut, label="Extrapolation")
    display(fig)
    
    append!(Ms,M)
    append!(errorInterpIns,maximum(errorInterpIn))
    append!(errorInterpOuts,maximum(errorInterpOut))
end

fig, ax = scatter(Ms,errorInterpIns,label="Interpolation")
ylims!(1e-12, 1e1)
ax.xlabel = "M"
ax.xscale = log10
ax.ylabel = "Relative error"
ax.yscale = log10
ax.title = "Interpolation 1D error"
scatter!(ax, Ms,errorInterpOuts;label="Extrapolation")
lines!(ax, Ms,(Ms/16).^-2; label=L"M^{-2}", color="red", linestyle=:dash)
axislegend(ax, position = :rt)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")


## Test integration

println("---------------------------------------------------------------------")
println("Relative error computing integrals")
println("")

Ms = Vector{Float64}()
errorInteg_with_log_binss = Vector{Float64}()
errorInteg_quad_with_log_binss = Vector{Float64}()
for j = 4:16
    M = 2^j
    kmin = 1e-3
    kmax = 1e+4

    Nk = wave_spectrum(kmin,kmax,M)
    kk = Nk.kk
    λ = Nk.λ
    @. Nk.nk = nk_test(kk)

    inteExact = 1.0
    inte = integrate_with_log_bins(Nk)
    errorInteg_with_log_bins = abs(inte - inteExact) / inteExact
    inte = integrate_quad_with_log_bins(Nk,1,-1)
    errorInteg_quad_with_log_bins = abs(inte - inteExact) / inteExact
    
    println("M = ", M, ": errors with_log_bins: ", errorInteg_with_log_bins, ", error quad_with_log_bins: ", errorInteg_quad_with_log_bins)

    append!(Ms,M)
    append!(errorInteg_with_log_binss,errorInteg_with_log_bins)
    append!(errorInteg_quad_with_log_binss,errorInteg_quad_with_log_bins)
end
fig, ax = scatter(Ms,errorInteg_with_log_binss,label="with_log_bins")
ylims!(1e-15, 1e1)
ax.xlabel = "M"
ax.xscale = log10
ax.ylabel = "Relative error"
ax.yscale = log10
ax.title = "Integration 1D error"
scatter!(ax,Ms,errorInteg_quad_with_log_binss;label="quad_with_log_bins")
lines!(ax, Ms,(Ms/16).^-2; label=L"M^{-2}", color="red", linestyle=:dash)
axislegend(ax, position = :rt)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")



println("---------------------------------------------------------------------")
println("Relative error computing integrals for linear grid")
println("")

Ms = Vector{Float64}()
errorIntegs = Vector{Float64}()
for j = 4:16
    M = 2^j 
    kmin = -1
    kmax = 1
    kk = vcat(LinRange(kmin,kmax,M))

    Nk = field_grid_1D(kk)
    kk = Nk.kk

    @. Nk.F = nk_test4(kk)

    inteExact = 1.0
    inte = integrate(Nk)
    errorInteg = abs(inte - inteExact) / inteExact
    
    println("M = ", M, ": errors trapz: ", errorInteg)

    append!(Ms,M)
    append!(errorIntegs,errorInteg)
end
fig, ax = scatter(Ms,errorIntegs,label="with_grid")
ylims!(1e-15, 1e1)
ax.xlabel = "M"
ax.xscale = log10
ax.ylabel = "Relative error"
ax.yscale = log10
ax.title = "Integration 1D error"
lines!(ax, Ms,(Ms/16).^-2; label=L"M^{-2}", color="red", linestyle=:dash)
axislegend(ax, position = :rt)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")


#benchmarking
#using BenchmarkTools
#@btime
