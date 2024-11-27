@doc raw"""
    plot_density!(Run; ρ=one, fig=nothing, ax=nothing, hm=nothing, label=L"$q_k$", title="Quantity spectrum", xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, xcomp = 0.)

Plot spectrum of a quantity whose spectral density is `ρ`. See [`density_spectrum!`](@ref).
* `fig`: Figure
* `ax`: Axis of the figure
* `hm`: Colormap of the figure
* `xcomp`: spectral slope compensation 

You can overwrite on a plot by specifying `fig`, `ax` and `hm`. Otherwise, it will create a new plot. Return `fig`, `ax` and `hm`.

"""
function plot_density!(Run; ρ=one, fig=nothing, ax=nothing, hm=nothing, label=L"$q_k$", title="Quantity spectrum", xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, xcomp = 0.)
  if Run.Nk_arguments == 1
    kk = Run.Nk.kk
    qk = Run.F1.nk
    density_spectrum!(qk, Run, ρ)
    fig, ax = plot_1D_base!(kk, abs.(qk) .* kk .^xcomp; fig=fig, ax=ax, xlabel=L"$k$", ylabel=label, title=title, xlims=xlims, ylims=ylims, xscale=xscale, yscale=yscale)
    hm = nothing
  elseif Run.Nk_arguments == 2
    kkh = Run.Nk.kkh
    kkz = Run.Nk.kkz
    qk = Run.F1.nk
    density_spectrum!(qk, Run, ρ)
    fig, ax, hm = plot_heatmap_base!(kkh, kkz, qk; fig=fig, ax=ax, hm=hm, xlabel=L"$k_h$", ylabel=L"$k_z$", title=title, xlims=xlims, ylims=ylims, zlims=zlims, xscale=xscale, yscale=yscale, zscale=zscale)
  else
    @error "wave action with three arguments not yet implemented"
  end
  return fig, ax, hm
end

@doc raw"""
    function plot_wave_action!(Run; fig=nothing, ax=nothing, hm=nothing, xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, xcomp = 0.)

Plot wave action spectrum. It is [`plot_density!`](@ref) with `ρ=1`.

"""
function plot_wave_action!(Run; fig=nothing, ax=nothing, hm=nothing, xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, xcomp = 0.)
  fig, ax, hm = plot_density!(Run; ρ=one, fig=fig, ax=ax, hm=hm, label=L"$n_k$", title="Wave action spectrum", xlims=xlims, ylims=ylims, zlims=zlims, xscale=xscale, yscale=yscale, zscale=zscale, xcomp=xcomp)
  return fig, ax, hm
end

@doc raw"""
    function plot_energy!(Run; fig=nothing, ax=nothing, hm=nothing, xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, xcomp = 0.)

Plot energy spectrum. It is [`plot_density!`](@ref) with `ρ=ω`.

"""
function plot_energy!(Run; fig=nothing, ax=nothing, hm=nothing, xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, xcomp = 0.)
  fig, ax, hm = plot_density!(Run; ρ=Run.ω, fig=fig, ax=ax, hm=hm, label=L"$e_k$", title="Energy spectrum", xlims=xlims, ylims=ylims, zlims=zlims, xscale=xscale, yscale=yscale, zscale=zscale, xcomp=xcomp)
  return fig, ax, hm
end

@doc raw"""
    function plot_1D_spectra_density!(Run; ρ=one, fig=nothing, ax=nothing, ylabel=L"$q_k$", title="Quantity 1D spectra", xlims=nothing, ylims=nothing, xscale=log10, yscale=log10, xcomp = 0.)

Plot ``k``, ``k_h`` and ``k_z`` spectra of a quantity whose spectral density is `ρ`. See [`isotropic_density_spectrum!`](@ref), [`kh_density_spectrum!`](@ref) and [`kz_density_spectrum!`](@ref).
* `fig`: Figure
* `ax`: Axis of the figure
* `xcomp`: spectral slope compensation 

You can overwrite on a plot by specifying `fig` and `ax`. Otherwise, it will create a new plot. Return `fig` and `ax`.

"""
function plot_1D_spectra_density!(Run; ρ=one, fig=nothing, ax=nothing, ylabel=L"$q_k$", title="Quantity 1D spectra", xlims=nothing, ylims=nothing, xscale=log10, yscale=log10, xcomp = 0.)
  if Run.Nk_arguments == 1
    @error "No kh and kz spectra for isotropic systems. Use density_spectrum(Run; ρ)."
  elseif Run.Nk_arguments == 2
    kkh = Run.Nk.kkh
    kkz = Run.Nk.kkz

    kk, q_k = isotropic_density_spectrum!(Run, ρ)
    q_kh = kh_density_spectrum!(Run, ρ)
    q_kz = kz_density_spectrum!(Run, ρ)

    fig, ax = plot_1D_base!(kk, q_k .* kk.^xcomp; fig=fig, ax=ax, xlabel=L"$k, k_h, k_z$", ylabel=ylabel, title=title, xlims=xlims, ylims=ylims, xscale=xscale, yscale=yscale)
    lines!(ax, kkh, q_kh .* kk.^xcomp, color = "red", label=L"vs $k_h$")
    lines!(ax, kkz, q_kz .* kk.^xcomp, color = "magenta", label=L"vs $k_z$")

    axislegend(position = :cb)
  else
    @error "wave action with three arguments not yet implemented"
  end
  return fig, ax
end

@doc raw"""
    plot_2D_slices_khkz!(Run; ρ=one, dir=1, nslices=8, fig=nothing, ax=nothing, cb=nothing, zlabel="", title=L"Spectrum slices", zlims=nothing, colormap=:inferno)

Plot ``k_h`` (`dir=1`) or ``k_z`` (`dir=2`) slices of the spectrum of the quantity with spectral density `ρ`. See [`density_spectrum!`](@ref) and [`plot_2D_slices_base!`](@ref).

"""
function plot_2D_slices_khkz!(Run; ρ=one, dir=1, nslices=8, fig=nothing, ax=nothing, cb=nothing, zlabel="", title=L"Spectrum slices", zlims=nothing, colormap=:inferno)
  kkh = Run.Nk.kkh
  kkz = Run.Nk.kkz
  qk = zeros(Run.Nk.Mh,Run.Nk.Mz)
  density_spectrum!(qk, Run, ρ)
  if dir == 1
    clabel = L"$k_h$"
    clims = (kkh[1], kkh[end])
    slabel = L"$k_z$"
    slims = (kkz[1], kkz[end])
  elseif dir == 2
    clabel = L"$k_z$"
    clims = (kkz[1], kkz[end])
    slabel = L"$k_h$"
    slims = (kkh[1], kkh[end])
  else
    @error "dir must be 1 (for kh slices) or 2 (for kz slices)"
  end
  fig, ax, cb = plot_2D_slices_base!(kkh, kkz, qk; dir=dir, nslices=nslices, fig=fig, ax=ax, cb=cb, slabel=slabel, zlabel=zlabel, clabel=clabel, title=title, slims=slims, zlims=zlims, clims=clims, colormap=colormap)  
  
end

@doc raw"""
    plot_2D_slices(Run; ρ=1, xlims=nothing, ylims=nothing, xscale=log10, yscale=log10, sthks=0:0.2:1)

Plot slices of the spectrum of the quantity with spectral density `ρ` at constant angles `sthks` × π/2 with the vertical axis ``{\bf e}_z``. Return the figure `fig` and its axis `ax`.

"""
function plot_2D_slices(Run; ρ=1, xlims=nothing, ylims=nothing, xscale=log10, yscale=log10, sthks=0:0.2:1)
  @assert Run.Nk_arguments == 2 
  if xlims === nothing
    xlims = (minimum(Run.Nk.kkh), maximum(Run.Nk.kkh))
  end

  if ylims === nothing
    ylims = (minimum(Run.Nk.nk)/1.1, maximum(Run.Nk.nk)*1.1)
  end
  
  if xscale == :log10
    xticks = LogRange(xlims[1],xlims[2],5)
  else
    xticks = LinRange(xlims[1],xlims[2],5)
  end

  if yscale == :log10
    yticks = LogRange(ylims[1],ylims[2],5)
  else
    yticks = LinRange(ylims[1],ylims[2],5)
  end

  interp = Run.interp_scheeme
  WavKinS.update_coeff_interp!(interp, Run.Nk)
  kk = Run.Nk.kkh
  sk = zeros(length(kk))

  fig = Figure(size = (400, 300))  
  ax = Axis(fig[1, 1],
    xlabel = L"k",
    xlabelsize=fontsizeAxisLabel,
    xscale=xscale,
    xticklabelsize=fontsizeAxis,
    ylabel=L"\rho n(k, \theta_{\bf k})",
    ylabelsize=fontsizeAxisLabel,
    yscale=yscale,
    yticklabelsize=fontsizeAxis,
    title="Spectrum slices",
    titlesize=fontsizeTitle,
  )
  
  for i in eachindex(sthks)
    kkh = kk * sin(pi*sthks[i]/2)
    kkz = kk * cos(pi*sthks[i]/2)
    for ik in eachindex(kk)
      sk[ik] = ρ(kkh[ik], kkz[ik]) * val_nk(Run.interp_scheeme, Run.Nk, kkh[ik], kkz[ik])
    end
    lines!(ax, kk, sk, colormap=:inferno, color=sthks[i], colorrange=(0,1.25), label=L"$2 \theta_k / \pi = %$(round(sthks[i],digits=2))$")  
  end
  axislegend(ax, position = :lb, labelsize=10, nbanks=2)
  xlims!(ax,xlims)
  ylims!(ax,ylims)
  return fig, ax
end