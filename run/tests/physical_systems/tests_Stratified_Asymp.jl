push!(LOAD_PATH, "../")
using WavKinS
using GLMakie
using LaTeXStrings
include("../../../src/physical_systems/Stratified/basics.jl")



function nk_test1(kx, ky)
    return exp(-kx - abs(ky))
end

function nk_test2(kx, ky)
    return abs(ky)*exp(-kx - abs(ky)) * kx^1.5/(1. + abs(ky))
end

function nk_test3(kx, ky)
    return abs(ky^2)*exp(-kx - abs(ky)) * kx^1.5/(1. + abs(ky))/118
end

function nk_test4(kx, ky)
    return exp.(- ( (log10.(abs.(kx))/2+0).^2  + (log10.(abs.(ky))/2+0).^2)/0.05)
end

function nk_test(kx,ky)
    return nk_test3(kx,ky)
end


#=
println("---------------------------------------------------------------------")
println("Testing functions and properties of the discriminant and interaction coefficients")
println("")

kh = 1.0
kz = 1.0

ϵ = 1e-6 # threeshold to avoid p=±kh and q=0 for which k1z or k2z is 0, then ω1 or ω2 and g' are not defined
p = LinRange(-kh+ϵ,kh-ϵ,256)
q = LinRange(ϵ,512*kh,8192)

P,Q = meshgrid(p,q)
KH = kh*ones(length(p),length(q))
KZ = kz*ones(length(p),length(q))

K1H = (KH+P+Q)/2
K2H = (KH-P+Q)/2

############################################
# Resonant manifold (k = k1 + k2) in (p,q) #
############################################
# Branch +
K1Zp = k1zkp_vs.(KH,KZ,K1H,K2H)
K2Z = KZ - K1Zp
dω = ω_Stratified_Asymp.(KH,KZ) - ω_Stratified_Asymp.(K1H,K1Zp) - ω_Stratified_Asymp.(K2H,K2Z)
@assert maximum(abs.(dω)) < 1e-9
# Branch -
K1Zm = k1zkm_vs.(KH,KZ,K1H,K2H)
K2Z = KZ - K1Zm
dω = ω_Stratified_Asymp.(KH,KZ) - ω_Stratified_Asymp.(K1H,K1Zm) - ω_Stratified_Asymp.(K2H,K2Z)
@assert maximum(abs.(dω)) < 1e-9


############################################
# Resonant manifold (k1 = k + k2) in (p,q) #
############################################
# Branch +
K1Zp = k1z1p_vs.(KH,KZ,K1H,K2H)
K2Z = K1Zp - KZ
dω = ω_Stratified_Asymp.(K1H,K1Zp) - ω_Stratified_Asymp.(KH,KZ) - ω_Stratified_Asymp.(K2H,K2Z)
@assert maximum(abs.(dω)) < 1e-9
# Branch -
K1Zm = k1z1m_vs.(KH,KZ,K1H,K2H)
K2Z = K1Zm - KZ
dω = ω_Stratified_Asymp.(K1H,K1Zm) - ω_Stratified_Asymp.(KH,KZ) - ω_Stratified_Asymp.(K2H,K2Z)
@assert maximum(abs.(dω)) < 1e-9



println("---------------------------------------------------------------------")
println("Testing grids")
println("")

Mh = 7
Mz = 8
khmin = 3e-1
khmax = 1.3e1
kzmin = 1e-1
kzmax = 1.9e1

Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
Run = Stratified_Asymp(Nk)

kkh = Nk.kkh
kkz = Nk.kkz
kk = Nk.kk

for ih in eachindex(kkh)
    aa = Run.kin_box.aa[ih]
    ppp = kkh[ih]*ones(length(aa)) - aa
    ppm = aa - kkh[ih]*ones(length(aa))
    pp = vcat(ppm,ppp)
    qq = Run.kin_box.qq
    Mp = length(pp[end])
    Mq = length(qq)
    qq = vcat([0],qq)
    PP, QQ = meshgrid(pp,qq)

    fig, ax = scatter(PP[:],QQ[:])
    lines!(ax, [-kkh[ih], -kkh[ih]], [0, 2*khmax])
    lines!(ax, [kkh[ih], kkh[ih]], [0, 2*khmax])
    lines!(ax, [-kkh[ih], kkh[ih]], [2*khmax, 2*khmax])
    lines!(ax, [-kkh[ih], kkh[ih]], [khmin/Mq, khmin/Mq])
    lines!(ax, [-kkh[ih], kkh[ih]], [0, 0])
    lines!(ax, [-kkh[ih]+khmin/Mh, -kkh[ih]+khmin/Mh], [0, 2*khmax])
    lines!(ax, [kkh[ih]-khmin/Mh, kkh[ih]-khmin/Mh], [0, 2*khmax])
    ax.xlabel = L"p"
    ax.xscale = identity
    ax.ylabel = L"q"
    ax.yscale= identity
    ax.title = "Kinematic box Stratified_Asymp"
    display(fig)
    sleep(5)
end

println("")
println("---------------------------------------------------------------------")


println("---------------------------------------------------------------------")
println("Testing collisional integral for simple integrand")
println("")

for j in 3:6
    Mh = 2^j
    Mz = 2^(j+1)
    khmin = 1e-3
    khmax = 1e3
    kzmin = 1e-3
    kzmax = 1e3
 
    Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    Run = Stratified_Asymp(Nk)

    qmin = Run.kin_box.qq[1]
    qmax = Run.kin_box.qq[end]

    kkh = Nk.kkh
    kkz = Nk.kkz
    kk = Nk.kk
    KKH = kkh .* ones(length(kkz))'
    KKZ = ones(length(kkh)) .* kkz'
    ONE = ones(length(kkh),length(kkz))
    @. Nk.nk = nk_test.(KKH, KKZ);

    WavKinS.St_k!(Run)

    Flux = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    @. Flux.nk = Run.Sk.nk * WavKinS.ω_Stratified_Asymp.(KKH, KKZ);
    integ = WavKinS.integrate_with_log_bins_khkz()
    Sum_Flux = integrate(integ, Flux)
    
    factor_Sk = 4*pi

    ###########################################################
    # Check integration on integral Stk with L = 1
    Sk_theo = factor_Sk * 4*pi*(asinh.(sqrt.(qmax./(2*KKH)))) 
    # Check integration on integral Stk with F = q
    #Sk_theo = factor_Sk * 2*pi*(-2*KKH.*asinh.(sqrt.(qmax./(2*KKH))) + sqrt.(qmax*(2*KKH + qmax*ONE)))  
    # Check integration on integral Stk with F = p^2
    #Sk_theo = factor_Sk * 2*pi*(KKH.^2).*(asinh.(sqrt.(qmax./(2*KKH))))  
    ###########################################################

    println("Mh = ", Mh, ", Mz = ", Mz , " rel. error ", maximum(abs.((Run.Sk.nk-Sk_theo)./Sk_theo)))

    fig, ax = plot_heatmap_base!(log10.(kkh), log10.(kkz), abs.((Run.Sk.nk-Sk_theo)./Sk_theo), 
        xlabel=L"k_h",
        ylabel=L"k_z",
        title=L"Relative error $St_{\mathbf{k}}$",
        xscale=identity,
        yscale=identity,
        zscale=identity
    )
    display(fig)
    sleep(5)
end

println("")
println("---------------------------------------------------------------------")
=#

println("---------------------------------------------------------------------")
println("Testing collisional integral for energy conservation")
println("")

Ms = Vector{Float64}()
errorFluxs = Vector{Float64}()
for j = 3:7
    Mh = 2^j
    Mz = 2^j
    khmin = 1e-2
    khmax = 1e2
    kzmin = 1e-2
    kzmax = 1e2
 
    Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    Run = Stratified_Asymp(Nk; interp_scheeme=WavKinS.bilin_interp_khkz, integ_scheeme=integrate_with_log_bins_khkz)

    kkh = Nk.kkh
    kkz = Nk.kkz
    kk = Nk.kk
    KKH = kkh .* ones(length(kkz))'
    KKZ = ones(length(kkh)) .* kkz'
    @. Nk.nk = nk_test.(KKH, KKZ);

    WavKinS.St_k!(Run)


    Flux = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    @. Flux.nk = Run.Sk.nk * Run.ω.(KKH, KKZ) * KKH * Run.dΩ;

    Ene = energy(Run)
    FluxNumH = total_density_flux(Run, Run.ω)
    
    
    println("Mh = ", Mh, ", Mz = ", Mz, ", Integral flux num: dH/H=", FluxNumH / Ene)


    fig, ax = plot_heatmap_base!(log10.(kkh),log10.(kkz),Flux.nk,
        xlabel=L"$k_h$",
        ylabel=L"$k_z$",
        title=L"$\omega_{\mathbf{k}} St_{\mathbf{k}}$",
        xscale=identity,
        yscale=identity,
        zscale=identity
    )
    display(fig)
    sleep(5)
    append!(Ms,Mh)
    append!(errorFluxs,abs(FluxNumH / Ene))
end
fig, ax = scatter(Ms,errorFluxs,label="with_log_bins")
ylims!(1e-4, 1e0)
ax.xlabel = L"M"
ax.xscale = log10
ax.ylabel = "Sum energy flux"
ax.yscale = log10
ax.title = "Energy flux Stratified_Asymp error"
lines!(ax, Ms, 1e2 * Ms.^-2; label=L"M^{-2}", color="red", linestyle=:dash)
axislegend(ax, position = :rt)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")
