# Tests for diagnostics

push!(LOAD_PATH, "../")
using WavKinS
using GLMakie


# Test isotropic_density_spectrum!
Mh = 128
Mz = 256
khmin = 1e-3
khmax = 1e0
kzmin = 1e-3
kzmax = 1e0
Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
Run = Acoustic2D_khkz(Nk);


KH, KZ = meshgrid(Run.Nk.kkh,Run.Nk.kkz)
K = sqrt.(KH.^2 + KZ.^2)
α =2
@. Run.Nk.nk = K^α

kk, nk_isotropic = isotropic_density_spectrum!(Run,Run.ω)
inds = nk_isotropic .>0
fig, ax = plot_1D_base!(kk[inds],nk_isotropic[inds]; xlabel=L"$k$", ylabel=L"$n(k)$", title="Test isotropic_density_spectrum", ylims=(1e-12,1e3))
plot_theo!(ax, kk, kk.^(α+2))

display(fig)
sleep(5)





# Test kh_density_spectrum! and kz_density_spectrum!
Mh = 1024
Mz = 1024
khmin = 1e-3
khmax = 1e0
kzmin = 1e-3
kzmax = 1e0
Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
Run = Stratified_Asymp(Nk);
kkh = Nk.kkh
kkz = Nk.kkz

KH, KZ = meshgrid(Run.Nk.kkh,Run.Nk.kkz)
αh = -0.5
αz = 2
@. Run.Nk.nk = KH.^αh .* KZ.^αz

spectrum_vs_kh = kh_density_spectrum!(Run,WavKinS.one)
fig, ax = plot_1D_base!(Nk.kkh,spectrum_vs_kh; xlabel=L"$k_h$", ylabel=L"$n(k_h)$", title="Test kh_density_spectrum", ylims=(1e-2,1e1))
factor = (kkz[end]^(αz+1) - kkz[1]^(αz+1)) / (αz+1)
plot_theo!(ax, kkh, kkh.^(αh + Run.dimension - 1) * Run.dΩ * factor) 
display(fig)
sleep(5)

spectrum_vs_kz = kz_density_spectrum!(Run,WavKinS.one)
fig, ax = plot_1D_base!(kkz,spectrum_vs_kz; xlabel=L"$k_z$", ylabel=L"$n(k_z)$", title="Test kz_density_spectrum", ylims=(1e-6,1e2))
factor = (kkh[end]^(αh + Run.dimension) - kkh[1]^(αh + Run.dimension)) / (αh + Run.dimension)
plot_theo!(ax, kkz, kkz.^αz * Run.dΩ * factor)
display(fig)
sleep(5)