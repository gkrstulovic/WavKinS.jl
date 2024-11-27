# Special outputs and diagnostics are defined here


function ξ_vs(kh,kz)
  return kh/kz^2  
end

@doc raw"""
  `compute_ωkz_spectrum(Run; ρ=one)`

Compute ``(\omega_{\bf k},|k_z|)`` spectrum for the quantity with spectral density `ρ`. It relates to the 3D wave action spectrum ``n_{\bf k}`` by ``s(\omega_{\bf k},|k_z|) = 2 \rho_{\bf k} \omega_{\bf k} (k_z/N)^2 \mathrm{d} \Omega n_{\rm k}``.
"""
function compute_ωkz_spectrum(Run; ρ=one)
  @assert Run.Nk_arguments == 2 

  kkh = Run.Nk.kkh
  kkz = Run.Nk.kkz
  interp = Run.interp_scheeme
  WavKinS.update_coeff_interp!(interp, Run.Nk)
  ωmin = Run.ω(kkh[1],kkz[end])
  ωmax = Run.ω(kkh[end],kkz[1])
  ωk = LogRange(ωmin,ωmax,Run.Nk.Mh)
  qk = zeros(length(ωk),Run.Nk.Mz)

  for iω in eachindex(ωk), iz in eachindex(kkz)
    kh = kkz[iz] * ωk[iω] / Run.N
    kz = kkz[iz]
    if kh >= kkh[1] && kh <= kkh[end] && kz >= kkz[1] && kz <= kkz[end] # in order to avoid plotting extrapolation errors
        qk[iω,iz] = 2 * ρ(kh,kz) * ωk[iω] * (kz/Run.N)^2 * Run.dΩ * val_nk(Run.interp_scheeme, Run.Nk, kh, kz) 
    else
        qk[iω,iz] = NaN
    end
  end

  return ωk, kkz, qk
end

@doc raw"""
  `compute_ξkz_spectrum(Run; ρ=one)`

Compute ``(ξ_{\bf k},|k_z|)`` spectrum for the quantity with spectral density `ρ`, where ``ξ_{\rm k}=k_h/k_z^2``. It relates to the 3D wave action spectrum ``n_{\bf k}`` by ``s(ξ_{\bf k},|k_z|) = 2 \rho_{\bf k} ξ_{\bf k} k_z^4 \mathrm{d} \Omega n_{\rm k}``.
"""
function compute_ξkz_spectrum(Run; ρ=one)
  @assert Run.Nk_arguments == 2 

  kkh = Run.Nk.kkh
  kkz = Run.Nk.kkz
  interp = Run.interp_scheeme
  WavKinS.update_coeff_interp!(interp, Run.Nk)
  ξmin = ξ_vs(kkh[1],kkz[end])
  ξmax = ξ_vs(kkh[end],kkz[1])
  ξk = LogRange(ξmin,ξmax,Run.Nk.Mh)
  qk = zeros(length(ξk),Run.Nk.Mz)

  for iξ in eachindex(ξk), iz in eachindex(kkz)
    kz = kkz[iz]
    kh = ξk[iξ] * kz^2
    if kh >= kkh[1] && kh <= kkh[end] && kz >= kkz[1] && kz <= kkz[end] # in order to avoid plotting extrapolation errors
        qk[iξ,iz] = 2* ρ(kh,kz) * ξk[iξ] * kz^4 * Run.dΩ * val_nk(Run.interp_scheeme, Run.Nk, kh, kz) 
    else
        qk[iξ,iz] = NaN
    end
  end

  return ξk, kkz, qk
end

@doc raw"""
  `plot_ωkz_spectrum!(Run; ρ=one, title="Quantity spectrum")`

Plot ``(ω_{\bf k},|k_z|)`` spectrum for the quantity with spectral density `ρ`. See [`compute_ωkz_spectrum`](@ref).
"""
function plot_ωkz_spectrum!(Run; ρ=one, fig=nothing, ax=nothing, hm=nothing, title="Quantity spectrum", zlims=nothing, xscale=log10, yscale=log10, zscale=log10, colormap=:jet)
  ωk, kkz, qk = compute_ωkz_spectrum(Run; ρ=ρ)
  fig, ax, hm = plot_heatmap_base!(ωk, kkz, qk; fig=fig, ax=ax, hm=hm, xlabel=L"$\omega_k$", ylabel=L"$k_z$", title=title, zlims=zlims, xscale=xscale, yscale=yscale, zscale=zscale, colormap=colormap)
  return fig, ax, hm
end

@doc raw"""
  `plot_ξkz_spectrum!(Run; ρ=one, title="Quantity spectrum")`

Plot ``(ξ_{\bf k},|k_z|)`` spectrum for the quantity with spectral density `ρ`. See [`compute_ξkz_spectrum`](@ref).
"""
function plot_ξkz_spectrum!(Run; ρ=one, fig=nothing, ax=nothing, hm=nothing, title="Quantity spectrum", zlims=nothing, xscale=log10, yscale=log10, zscale=log10, colormap=:jet)
  ξk, kkz, qk = compute_ξkz_spectrum(Run; ρ=ρ)
  fig, ax, hm = plot_heatmap_base!(ξk, kkz, qk; fig=fig, ax=ax, hm=hm, xlabel=L"$ξ_k$", ylabel=L"$k_z$", title=title, zlims=zlims, xscale=xscale, yscale=yscale, zscale=zscale, colormap=colormap)
  return fig, ax, hm
end

@doc raw"""
  `plot_slices_ωkz_spectrum!(Run; ρ=one, title="Quantity spectrum")`

Plot slices of the ``(ω_{\bf k},|k_z|)`` spectrum for the quantity with spectral density `ρ`. See [`compute_ωkz_spectrum`](@ref).
"""
function plot_slices_ωkz_spectrum!(Run; ρ=one, nslices=8, fig=nothing, ax=nothing, cb=nothing, zlabel=L"n(\omega,|k_z|)", title=L"Spectrum slices", kzlims=nothing, zlims=nothing)
  ωk, kkz, qk = compute_ωkz_spectrum(Run; ρ=ρ)
  fig, ax, cb = plot_2D_slices_base!(ωk, kkz, qk; dir=1, nslices=nslices, fig=fig, ax=ax, cb=cb, slabel=L"$k_z$", zlabel=zlabel, clabel=L"$\omega_k$", title=title, slims=kzlims, zlims=zlims, sscale=log10, zscale=log10, cscale=log10)
  return fig, ax, cb
end

@doc raw"""
  `plot_slices_kzω_spectrum!(Run; ρ=one, title="Quantity spectrum")`

Plot slices of the ``(ω_{\bf k},|k_z|)`` spectrum for the quantity with spectral density `ρ`. See [`compute_ωkz_spectrum`](@ref).
"""
function plot_slices_kzω_spectrum!(Run; ρ=one, nslices=8, fig=nothing, ax=nothing, cb=nothing, zlabel=L"n(\omega,|k_z|)", title=L"Spectrum slices", zlims=nothing, ωlims=nothing)
  ωk, kkz, qk = compute_ωkz_spectrum(Run; ρ=ρ)
  fig, ax, cb = plot_2D_slices_base!(ωk, kkz, qk; dir=2, nslices=nslices, fig=fig, ax=ax, cb=cb, slabel=L"$\omega_k$", zlabel=zlabel, clabel=L"$k_z$", title=title, slims=ωlims, zlims=zlims, sscale=log10, zscale=log10, cscale=log10)
  return fig, ax, cb
end

@doc raw"""
  `plot_slices_ξkz_spectrum!(Run; ρ=one, title="Quantity spectrum")`

Plot slices of the ``(ξ_{\bf k},|k_z|)`` spectrum for the quantity with spectral density `ρ`. See [`compute_ξkz_spectrum`](@ref).
"""
function plot_slices_ξkz_spectrum!(Run; ρ=one, nslices=32, fig=nothing, ax=nothing, cb=nothing, zlabel=L"n(\xi,|k_z|)", title=L"Spectrum slices", kzlims=nothing, zlims=nothing)
  ξk, kkz, qk = compute_ξkz_spectrum(Run; ρ=ρ)
  fig, ax, cb = plot_2D_slices_base!(ξk, kkz, qk; dir=1, nslices=nslices, fig=fig, ax=ax, cb=cb, slabel=L"$k_z$", zlabel=zlabel, clabel=L"$\xi_k$", title=title, slims=kzlims, zlims=zlims, sscale=log10, zscale=log10, cscale=log10)
  return fig, ax, cb
end

@doc raw"""
  `plot_slices_kzξ_spectrum!(Run; ρ=one, title="Quantity spectrum")`

Plot slices of the ``(ξ_{\bf k},|k_z|)`` spectrum for the quantity with spectral density `ρ`. See [`compute_ξkz_spectrum`](@ref).
"""
function plot_slices_kzξ_spectrum!(Run; ρ=one, nslices=32, fig=nothing, ax=nothing, cb=nothing, zlabel=L"n(\xi,|k_z|)", title=L"Spectrum slices", zlims=nothing, ξlims=nothing)
  ξk, kkz, qk = compute_ξkz_spectrum(Run; ρ=ρ)
  fig, ax, cb = plot_2D_slices_base!(ξk, kkz, qk; dir=2, nslices=nslices, fig=fig, ax=ax, cb=cb, slabel=L"$\xi_k$", zlabel=zlabel, clabel=L"$k_z$", title=title, slims=ξlims, zlims=zlims, sscale=log10, zscale=log10, cscale=log10)
  return fig, ax, cb
end

@doc raw"""
  `plot_ω_spectrum_density!(Run; ρ=one, title="Quantity spectrum")`

Plot frequency spectrum for the quantity with spectral density `ρ` by integrating the ``(ω_{\bf k},|k_z|)`` spectrum over ``|k_z|``. See [`compute_ωkz_spectrum`](@ref). 
"""
function plot_ω_spectrum_density!(Run; ρ=one, fig=nothing, ax=nothing, ylabel=L"n(\omega_{\bf k})", title=L"Frequency spectrum", ωlims=nothing, ylims=nothing, xcomp=0.0)
  ωk, kkz, qkk = compute_ωkz_spectrum(Run; ρ=ρ)
  q_ω = zeros(length(ωk))
  inds = .!(qkk .!== NaN)
  @. qkk[inds] = 0.0
  for iω in eachindex(ωk)
    q_ω[iω] = integrate_with_log_bins(Run.Nk.Mz, Run.Nk.λz, Run.Nk.logλz, kkz, qkk[iω,:])
  end
  fig, ax = plot_1D_base!(ωk, q_ω .* ωk.^xcomp; fig=fig, ax=ax, xlabel=L"$\omega_{\bf k}$", ylabel=ylabel, title=title, xlims=ωlims, ylims=ylims, xscale=log10, yscale=log10)
  return fig, ax
end