@doc raw"""
    plot_2D_flux(Nk::wave_spectrum_khkz, Fh::Matrix{Float64}, Fz::Matrix{Float64}; xlabel=L"$log_{10} k_h$", ylabel=L"$log_{10} k_z$", title="Fluxes")
 
Plot fluxes for anisotropic systems.
* `Nk`: structure containing the grid
* `Fh`: flux in horizontal wave vectors 
* `Fz`: flux in vertical wave vectors

"""
function plot_2D_flux(Nk::wave_spectrum_khkz, Fh::Matrix{Float64}, Fz::Matrix{Float64}; xlabel=L"$log_{10} k_h$", ylabel=L"$log_{10} k_z$", title="Fluxes")
  kkh = Nk.kkh
  kkz = Nk.kkz

  FFh = wave_spectrum_khkz(kkh[1],kkh[end],Nk.Mh,kkz[1],kkz[end],Nk.Mz)
  @. FFh.nk = Fh
  interph = bilin_interp_khkz(kkh, kkz)
  #interph = cpow_interp_khkz(kkh, kkz)
  update_coeff_interp!(interph, FFh)

  FFz = wave_spectrum_khkz(kkh[1],kkh[end],Nk.Mh,kkz[1],kkz[end],Nk.Mz)
  @. FFz.nk = Fz
  interpz = bilin_interp_khkz(kkh, kkz)
  #interpz = cpow_interp_khkz(kkh, kkz)
  update_coeff_interp!(interpz, FFz)

  
  f(x) = Point2f(
      val_nk(interph, FFh, exp10(x[1]), exp10(x[2])),
      val_nk(interpz, FFz, exp10(x[1]), exp10(x[2]))
  )

  fig, ax, pl = streamplot(f, (log10.(kkh[1]),log10.(kkh[end])), (log10.(kkz[1]),log10.(kkz[end])), 
    colormap = :jet, 
    gridsize=(32,32,32),
    colorscale=identity,
  )

  #fig, ax = arrows((1:10:Nk.Mh)/Nk.Mh, (1:10:Nk.Mz)./Nk.Mz, Fh[1:10:end,1:10:end],Fz[1:10:end,1:10:end])

  #Colorbar(fig[1, 2], pl)
  resize!(fig.scene,(400,300))
  ax.xlabel=xlabel
  ax.xlabelsize=fontsizeAxisLabel
  ax.xticklabelsize=fontsizeAxis
  ax.ylabel=ylabel
  ax.ylabelsize=fontsizeAxisLabel
  ax.yticklabelsize=fontsizeAxis
  ax.title=title
  ax.titlesize=fontsizeTitle

  xlims!(ax,log10(kkh[1]),log10(kkh[end]))
  ylims!(ax,log10(kkz[1]),log10(kkz[end]))
  return fig, ax
end

@doc raw"""
    plot_density_flux!(Run; ρ=one, fig=nothing, ax=nothing, label=L"$P_k$", title="Density flux", Disp = 1., xlims=nothing, ylims=nothing)

Plot fluxes of a quantity whose spectral density is `ρ`. See [`density_flux!`](@ref) and [`plot_2D_flux`](@ref).
* `fig`: Figure
* `ax`: Axis of the figure

You can overwrite on a plot by specifying `fig` and `ax`. Otherwise, it will create a new plot. Return `fig` and `ax`.

"""
function plot_density_flux!(Run; ρ=one, fig=nothing, ax=nothing, label=L"$P_k$", title="Density flux", Disp = 1., xlims=nothing, ylims=nothing)
  if Run.Nk_arguments == 1
    kk = Run.Nk.kk
    Pk = zeros(Run.Nk.M)
    density_flux!(Pk, Run, ρ)
    fig, ax = plot_1D_base!(kk, Pk./Disp; fig=fig, ax=ax, xlabel=L"$k$", ylabel=label, title=title, xlims=xlims, ylims=ylims, yscale=identity)
  elseif Run.Nk_arguments == 2
    Pkh = zeros(Run.Nk.Mh, Run.Nk.Mz)
    Pkz = zeros(Run.Nk.Mh, Run.Nk.Mz)
    density_flux!([Pkh, Pkz], Run, ρ)
    fig, ax = plot_2D_flux(Run.Nk, Pkh, Pkz; title=title)
  else
    @error "wave action with three arguments not yet implemented"
  end
  return fig, ax
end

@doc raw"""
    plot_wave_action_flux!(Run; fig=nothing, ax=nothing, Disp = 1., xlims=nothing, ylims=nothing)

Plot wave action fluxes. It is [`plot_density_flux!`](@ref) with `ρ=1`.

"""
function plot_wave_action_flux!(Run; fig=nothing, ax=nothing, Disp = 1., xlims=nothing, ylims=nothing)
  fig, ax = plot_density_flux!(Run; ρ=one, fig=fig, ax=ax, label=L"$P_k$", title="Wave action flux", Disp=Disp, xlims=xlims, ylims=ylims)
  return fig, ax
end

@doc raw"""
    plot_energy_flux!(Run; fig=nothing, ax=nothing, Disp = 1., xlims=nothing, ylims=nothing)

Plot energy fluxes. It is [`plot_density_flux!`](@ref) with `ρ=ω`.

"""
function plot_energy_flux!(Run; fig=nothing, ax=nothing, Disp = 1., xlims=nothing, ylims=nothing)
  fig, ax = plot_density_flux!(Run; ρ=Run.ω, fig=fig, ax=ax, label=L"$P_k$", title="Energy flux", Disp=Disp, xlims=xlims, ylims=ylims)
  return fig, ax
end

@doc raw"""
    plot_density_flux_isotropic!(Run; ρ=one, fig=nothing, ax=nothing, label=L"$P_k$", title="Isotropic quantity flux", Disp=1., xlims=nothing, ylims=nothing, xscale=log10, yscale=identity)

Plot the isotropic (after sum over angle) fluxes of a quantity whose spectral density is `ρ`. Corresponds to the standard flux for isotropic systems. See [`density_flux_isotropic!`](@ref).
* `fig`: Figure
* `ax`: Axis of the figure

You can overwrite on a plot by specifying `fig` and `ax`. Otherwise, it will create a new plot. Return `fig` and `ax`.

"""
function plot_density_flux_isotropic!(Run; ρ=one, fig=nothing, ax=nothing, label=L"$P_k$", title="Isotropic quantity flux", Disp=1., xlims=nothing, ylims=nothing, xscale=log10, yscale=identity)
  if Run.Nk_arguments == 1
    # We use the standard function for isotropic systems
    fig, ax = plot_density_flux!(Run; ρ=ρ, fig=fig, ax=ax, label=label, title=title, Disp = 1., xlims=xlims, ylims=ylims)
  elseif Run.Nk_arguments == 2
    kk = Run.Nk.kkh
    Pk = zeros(Run.Nk.Mh)
    density_flux_isotropic!(Pk, Run, ρ)
    fig, ax = plot_1D_base!(kk, Pk./Disp; fig=ax, ax=ax, xlabel=L"$k$", ylabel=label, title=title, xlims=xlims, ylims=ylims, xscale=xscale, yscale=yscale)
  else
    @error "wave action with three arguments not yet implemented"
  end
  return fig, ax
end

@doc raw"""
    plot_wave_action_flux_isotropic!(Run; fig=nothing, ax=nothing, Disp=1., xlims=nothing, ylims=nothing, xscale=log10, yscale=identity)

Plot the isotropic (after sum over angle) wave action fluxes. It is [`plot_density_flux_isotropic!`](@ref) with `ρ=1`.

"""
function plot_wave_action_flux_isotropic!(Run; fig=nothing, ax=nothing, Disp=1., xlims=nothing, ylims=nothing, xscale=log10, yscale=identity)
  fig, ax = plot_density_flux_isotropic!(Run; ρ=one, fig=fig, ax=ax, label=L"$P_k$", title="Isotropic wave action flux", Disp=Disp, xlims=xlims, ylims=ylims, xscale=xscale, yscale=yscale)
  return fig, ax
end

@doc raw"""
    plot_energy_flux_isotropic!(Run; fig=nothing, ax=nothing, Disp=1., xlims=nothing, ylims=nothing, xscale=log10, yscale=identity)

Plot the isotropic (after sum over angle) energy fluxes. It is [`plot_density_flux_isotropic!`](@ref) with `ρ=ω`.

"""
function plot_energy_flux_isotropic!(Run; fig=nothing, ax=nothing, Disp=1., xlims=nothing, ylims=nothing, xscale=log10, yscale=identity)
  fig, ax = plot_density_flux_isotropic!(Run; ρ=Run.ω, fig=fig, ax=ax, label=L"$P_k$", title="Isotropic energy flux", Disp=Disp, xlims=xlims, ylims=ylims, xscale=xscale, yscale=yscale)
  return fig, ax
end

@doc raw"""
    plot_density_flux_angular!(Run; ρ=one, fig=nothing, ax=nothing, label=L"$P_k$", title="Angular quantity flux", Disp = 1., xlims=nothing, ylims=nothing, xscale=identity, yscale=identity)

Plot the angular fluxes of a quantity whose spectral density is `ρ`. See [`density_flux_angular!`](@ref).
* `fig`: Figure
* `ax`: Axis of the figure

You can overwrite on a plot by specifying `fig` and `ax`. Otherwise, it will create a new plot. Return `fig` and `ax`.

"""
function plot_density_flux_angular!(Run; ρ=one, fig=nothing, ax=nothing, label=L"$P_k$", title="Angular quantity flux", Disp = 1., xlims=nothing, ylims=nothing, xscale=identity, yscale=identity)
  if Run.Nk_arguments == 1
    @error "No angular flux for isotropic systems. No plot."
  elseif Run.Nk_arguments == 2
    Mthk = Int(round(Run.Nk.Mh/2))
    thk = vec(0.0:pi/(2*(Mthk-1)):pi/2)
    Pthk = zeros(Mthk)
    density_flux_angular!(thk,Pthk,Run,ρ)
    fig, ax = plot_1D_base!(thk, Pthk./Disp; fig=fig, ax=ax, xlabel=L"$\theta_k$", ylabel=label, title=title, xlims=xlims, ylims=ylims, xscale=xscale, yscale=yscale)
  else
    @error "wave action with three arguments not yet implemented"
  end
  return fig, ax
end

@doc raw"""
    plot_wave_action_flux_angular!(Run; fig=nothing, ax=nothing, Disp = 1., xlims=nothing, ylims=nothing, xscale=identity, yscale=identity)

Plot the angular fluxes of wave action. It is [`plot_density_flux_angular!`](@ref) with `ρ=1`.

"""
function plot_wave_action_flux_angular!(Run; fig=nothing, ax=nothing, Disp = 1., xlims=nothing, ylims=nothing, xscale=identity, yscale=identity)
  fig, ax = plot_density_flux_angular!(Run; ρ=one, fig=fig, ax=ax, label=L"$P_k$", title="Angular wave action flux", Disp=Disp, xlims=xlims, ylims=ylims, xscale=xscale, yscale=yscale)
  return fig, ax
end

@doc raw"""
    plot_energy_flux_angular!(Run; fig=nothing, ax=nothing, Disp = 1., xlims=nothing, ylims=nothing, xscale=identity, yscale=identity)

Plot the angular fluxes of energy. It is [`plot_density_flux_angular!`](@ref) with `ρ=ω`.

"""
function plot_energy_flux_angular!(Run; fig=nothing, ax=nothing, Disp = 1., xlims=nothing, ylims=nothing, xscale=identity, yscale=identity)
  fig, ax = plot_density_flux_angular!(Run; ρ=Run.ω, fig=fig, ax=ax, label=L"$P_k$", title="Angular energy flux", Disp=Disp, xlims=xlims, ylims=ylims, xscale=xscale, yscale=yscale)
  return fig, ax
end
