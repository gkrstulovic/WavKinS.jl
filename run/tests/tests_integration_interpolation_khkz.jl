# Tests for integration and interpolation for anisotropic systems

push!(LOAD_PATH, "../")
using WavKinS
using GLMakie
using Random
using LaTeXStrings


function nk_test(kkh, kkz)
    # these functions have integrals equal to one (useful to test integration)
    nk = (kkh^2 * exp(-kkh / 50.0) / 250000) * (kkz^2 * exp(-kkz / 50.0) / 250000)
    return nk
end

function nk_test2(kkh, kkz)
    nk = kkh^(-1) * kkz^(-2) * exp(-((kkh^2 + kkz^2)^0.5/500.)^2)
    return nk
end

function nk_test3(kkh, kkz)
    nk = (kkh^2) * (kkz^2)
    return nk
end

function nk_test4(kkh, kkz)
    nk = exp( 1.0 +  2.5 * log(kkh)  + 0.5 *log(kkz) + 1.36 * log(kkh) * log(kkz) ) / (1.0 + exp( abs(kkh^2 +kkz^2) / 100)) 
    return nk
end

function nk_test5(kkh, kkz)
    nk = 0.5 + 1.0*kkh + 0.5*kkz + 3.365*kkh*kkz
    return nk
end

function n_test(kkh, kkz)
    return nk_test(kkh, kkz)
end




##   Test interpolation 

println("---------------------------------------------------------------------")
println("Relative error computing khkz interpolations")
println("")

Ms = Vector{Float64}()
errorInterpIns = Vector{Float64}()
errorInterpOuts = Vector{Float64}()
for j = 4:12
    Mh = 2^j
    Mz = 2^j
    khmin = 1e-2
    khmax = 1e+3
    kzmin = 1e-2
    kzmax = 1e+3

    Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    kkh = Nk.kkh
    kkz = Nk.kkz
    λh = Nk.λh
    λz = Nk.λz
    KKH = kkh .* ones(length(kkz))'
    KKZ = ones(length(kkh)) .* kkz'
    @. Nk.nk = n_test.(KKH, KKZ)

    interp = WavKinS.bilin_interp_khkz(kkh, kkz)
    #interp = WavKinS.cpow_interp_khkz(kkh, kkz)
    #interp = WavKinS.bilinlog_interp_khkz(kkh, kkz)
    interp = WavKinS.exp_interp_khkz(kkh, kkz)
    WavKinS.update_coeff_interp!(interp, Nk)

    InterpolationPoints = 500 * 2
    Random.seed!(1234)
    Ihrand = rand(InterpolationPoints) .* (Mh-1) .+ 1
    khtestIn = kkh[1] * λh .^ (Ihrand .- 1)
    Izrand = rand(InterpolationPoints) .* (Mz-1) .+ 1
    kztestIn = kkz[1] * λz .^ (Izrand .- 1)
    khtestOut = [khmin/2]
    kztestOut = [kzmin/2]
    NtestIn = similar(khtestIn)
    NtestOut = similar(khtestOut)

    for i in eachindex(khtestIn)
        NtestIn[i] = val_nk(interp, Nk, khtestIn[i], kztestIn[i])
    end
    for i in eachindex(khtestOut)
        NtestOut[i] = val_nk(interp, Nk, khtestOut[i], kztestOut[i])
    end

    inds = NtestIn .> 1e-40
    errorInterpIn = abs.((NtestIn[inds] .- n_test.(khtestIn[inds], kztestIn[inds])) ./ n_test.(khtestIn[inds], kztestIn[inds]))
    errorInterpOut = abs.((NtestOut .- n_test.(khtestOut, kztestOut)) ./ n_test.(khtestOut, kztestOut))

    println("Mh = ", Mh, ", Mz = ", Mz, ": relative error interpolation: ", maximum(errorInterpIn) , ", relative error extrapolation: ", maximum(errorInterpOut))
    fig1, ax1 = plot_surface_base!(Nk.kkh, Nk.kkz, Nk.nk; xlims=(0.4 * khmin, 1.1 * khmax), ylims=(0.4 * kzmin, 1.1 * kzmax), zlims=(1e-6, 1e7))
    scatter!(ax1, log10.(khtestIn[inds]), log10.(kztestIn[inds]), log10.(NtestIn[inds]), label="testIn")
    scatter!(ax1, log10.(khtestOut), log10.(kztestOut), log10.(NtestOut), label="testOut")
    #fig2, ax2 = scatter(log10.(khtestIn[inds]), log10.(errorInterpIn), marker_z=log10.(kztestIn[inds]), label="")
    #scatter!(ax2, log10.(khtestOut), log10.(errorInterpOut), marker_z=log10.(kztestOut), label="")
    sleep(5)
    display(fig1)
    append!(Ms,Mh)
    append!(errorInterpIns,maximum(errorInterpIn))
    append!(errorInterpOuts,maximum(errorInterpOut))
end
fig, ax = scatter(Ms,errorInterpIns,label="Interpolation")
ax.xlabel = "M"
ax.xscale = log10
ax.ylabel = "Relative error"
ax.yscale= log10
ax.title = "Interpolation 2D error"
scatter!(ax, Ms,errorInterpOuts;label="Extrapolation")
lines!(ax, Ms,(Ms/16).^-2; label=L"M^{-2}", color="red", linestyle=:dash)
axislegend(ax, position = :rt)
display(fig)
sleep(15)


## Test integration

println("---------------------------------------------------------------------")
println("Relative error computing khkz integrals")
println("")

Ms = Vector{Float64}()
errorInteg_with_log_binss = Vector{Float64}()
errorInteg_with_cpows = Vector{Float64}()
for j = 4:12
    Mh = 2^j
    Mz = 2^(j+1)
    khmin = 1e-2
    khmax = 1e+3
    kzmin = 1e-2
    kzmax = 1e+3
   
    Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    kkh = Nk.kkh
    kkz = Nk.kkz
    KKH = kkh .* ones(length(kkz))'
    KKZ = ones(length(kkh)) .* kkz'
    @. Nk.nk = n_test.(KKH, KKZ)
    
    #WavKinS.clean_waveaction!(Nk)

    inteExact = 1.0
    integ = WavKinS.integrate_with_log_bins_khkz()
    inte = integrate(integ, Nk, 1, -1, 1, -1)
    errorInteg_with_log_bins = abs(inte - inteExact) / inteExact
    integ = WavKinS.integrate_with_cpow_khkz() 
    inte = integrate(integ, Nk, 1, -1, 1, -1)
    errorInteg_with_cpow = abs(inte - inteExact) / inteExact

    println("Mh = ", Mh, ", Mz = ", Mz, ": errors integration with_log_bins: ", errorInteg_with_log_bins, ": errors integration with_cpow: ", errorInteg_with_cpow)

    append!(Ms,Mh)
    append!(errorInteg_with_log_binss, errorInteg_with_log_bins)
    append!(errorInteg_with_cpows, errorInteg_with_cpow)
end
fig, ax = scatter(Ms,errorInteg_with_log_binss,label="with_log_bins")
ax.xlabel = "M"
ax.xscale = log10
ax.ylabel = "Relative error"
ax.yscale = log10
ax.title = "Integration 2D error"
scatter!(ax,Ms,errorInteg_with_cpows;label="with_cpow")
lines!(ax, Ms,(Ms/16).^-2; label=L"M^{-2}", color="red", linestyle=:dash)
axislegend(ax, position = :rt)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")
