# Tests for plot

push!(LOAD_PATH, "../")
using WavKinS
using GLMakie


#=
# Check plot update
x = Vector(LinRange(0,2*pi,128))
Z = 0 * x .* x'
Nt = 256

fig, ax, hm = plot_surface_base!(x,x,Z; zlims=[-1, 1], xscale=identity, yscale=identity, zscale=identity, colormap=:deep)
display(fig)
for it in 0:Nt 
    t = it * 2*pi /Nt
    @. Z = sin(t) .* sin.(x .- t) .* sin.(x')
    plot_surface_base!(x,x,Z; fig=fig, ax=ax, hm=hm)
    sleep(0.1)
end
=#
#=
# Check plot_2D_slices_base
x = Vector(LinRange(0.1,2*pi,128))
y = Vector(LogRange(0.1,10000,128))
Z = (x.^-3) .* (y.^2)'
Nt = 256

fig, ax, cb = plot_2D_slices_base!(x,y,Z;dir=2,nslices=16)
display(fig)
for it in 0:Nt 
    t = it * 2*pi /Nt
    @. Z = Z = (x.^-cos(t)) .* (y.^(2*sin(t)))'
    plot_2D_slices_base!(x,y,Z; dir=2, nslices=16, fig=fig, ax=ax, cb=cb)
    sleep(0.1)
end
=#


# Check plot_energy_flux_isotropic
Mh = 1024
Mz = 512
khmin = 1e-2
khmax = 1e1
kzmin = 1e-2
kzmax = 1e1
Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
Run = Acoustic2D_khkz(Nk; a=0.0);

KH, KZ = meshgrid(Run.Nk.kkh,Run.Nk.kkz)
K = sqrt.(KH.^2 + KZ.^2)

α = -2
@. Run.Sk.nk = - K.^(α) 
@. Run.Sk.nk[K.>khmax] = 0.0

fig, ax = plot_energy_flux_isotropic!(Run)
kk = Nk.kkh
if α > -3.0
    plot_theo!(ax, kk, pi * kk.^(α+3) / (2*(α+3)))
elseif α == -3.0
    plot_theo!(ax, kk, pi * log.(sqrt(2)*kk/kk[1]) / 2 )
else
    println("Not integrable at zero")
end
display(fig)
sleep(5)

# Check plot_energy_flux_angular
Mh = 1024
Mz = 512
khmin = 1e-6
khmax = 1e0
kzmin = 1e-6
kzmax = 1e0
Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
Run = Acoustic2D_khkz(Nk; a=0.0);

KH, KZ = meshgrid(Run.Nk.kkh,Run.Nk.kkz)
K = sqrt.(KH.^2 + KZ.^2)
@. Run.Sk.nk =  -KH ./ K.^2
@. Run.Sk.nk[K.>khmax] = 0.0

fig, ax = plot_energy_flux_angular!(Run)
thk = vcat(0.0:pi/100:pi/2)
plot_theo!(ax, thk, sin.(thk)/2)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")


println("---------------------------------------------------------------------")
println("Test Anisotropic fluxes")
println("")

Mh = 128
Mz = 256
khmin = 1e-4
khmax = 1e0
kzmin = 1e-4
kzmax = 1e0
Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
Run = Stratified_Asymp(Nk);

KH, KZ = meshgrid(Run.Nk.kkh,Run.Nk.kkz)
K = sqrt.(KH.^2 + KZ.^2)
@. Run.Sk.nk =  -exp.(- ((KH .- 0.01).^2 + (KZ .- 0.01).^2) / 1e-6)
#@. Run.Sk.nk =  -KH.^(1)
#@. Run.Sk.nk[K.<=0.0] = NaN
#@. Run.Sk.nk[KH.<1e-5] = 1.0

fig, ax = plot_energy_flux!(Run)
#fig ,ax = heatmap(log10.(Nk.kkh), log10.(Nk.kkz), Run.Sk.nk)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")

