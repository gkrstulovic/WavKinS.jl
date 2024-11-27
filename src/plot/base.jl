using CairoMakie
using LaTeXStrings
using Printf

# Standard fontsizes
fontsizeAxis=12
fontsizeAxisLabel=18
fontsizeTitle=18


@doc raw"""
    plot_1D_base!(x, y; fig=nothing, ax=nothing, xlabel=L"$x$", ylabel=L"$y$", title="Title", xlims=nothing, ylims=nothing, xscale=log10, yscale=log10)

Basic plot of `y` vs `x`.
* `fig`: Figure
* `ax`: Axis of the figure

You can overwrite on a plot by specifying `fig` and `ax`. Otherwise, it will create a new plot. Return `fig` and `ax`.

"""
function plot_1D_base!(x, y; fig=nothing, ax=nothing, xlabel=L"$x$", ylabel=L"$y$", title="Title", xlims=nothing, ylims=nothing, xscale=log10, yscale=log10)  
  if fig === nothing || ax === nothing
    if xlims === nothing
      xlims = (minimum(x), maximum(x))
    end
  
    if ylims === nothing
      ylims = (minimum(y)/1.1, maximum(y)*1.1)
    end

    fig = Figure(size = (400, 300)) 
    ax = Axis(fig[1, 1],
      xlabel=xlabel,
      xlabelsize=fontsizeAxisLabel,
      xscale=xscale,
      xticklabelsize=fontsizeAxis,
      ylabel=ylabel,
      ylabelsize=fontsizeAxisLabel,
      yscale=yscale,
      yticklabelsize=fontsizeAxis,
      title=title,
      titlesize=fontsizeTitle,
    )

    xlims!(ax,xlims)
    ylims!(ax,ylims)
  end
  empty!(ax)
  lines!(ax,x,y,color="blue")
  return fig, ax
end

@doc raw"""
    plot_heatmap_base!(x, y, Z; fig=nothing, ax=nothing, hm=nothing, xlabel=L"$x$", ylabel=L"$y$", title="Title", xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, colormap=:jet)  

Basic heatmap plot of `Z` vs `x` and `y`.
* `fig`: Figure
* `ax`: Axis of the figure
* `hm`: Colormap of the figure

You can overwrite on a plot by specifying `fig`, `ax` and `hm`. Otherwise, it will create a new plot. Return `fig`, `ax` and `hm`.

"""
function plot_heatmap_base!(x, y, Z; fig=nothing, ax=nothing, hm=nothing, xlabel=L"$x$", ylabel=L"$y$", title="Title", xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, colormap=:jet)  
  if fig === nothing || ax === nothing || hm === nothing
    if xlims === nothing
      xlims = (minimum(x), maximum(x))
    end

    if ylims === nothing
      ylims = (minimum(y), maximum(y))
    end

    if zlims === nothing
      if zscale == log10
        inds = Z .> 0
        zlims = (minimum(Z[inds])/1.1, maximum(Z[inds])*1.1)
      else
        zlims = (minimum(Z)/1.1, maximum(Z)*1.1)
      end
    end

    fig, ax, hm = heatmap(x, y, Z; 
      colormap = colormap,
      colorrange = (zlims),
      colorscale = zscale, 
      axis = (; xscale = xscale, yscale = yscale),
    )
    Colorbar(fig[1, 2], hm)

    resize!(fig.scene,(400,300))
    ax.xlabel=xlabel
    ax.xlabelsize=fontsizeAxisLabel
    ax.xticklabelsize=fontsizeAxis
    ax.ylabel=ylabel
    ax.ylabelsize=fontsizeAxisLabel
    ax.yticklabelsize=fontsizeAxis
    ax.title=title
    ax.titlesize=fontsizeTitle

    xlims!(ax,xlims)
    ylims!(ax,ylims)
  else
    heatmap!(ax, x, y, Z; colormap=hm.colormap, colorrange=hm.colorrange, colorscale=hm.colorscale, overdraw=true)
  end
  return fig, ax, hm
end

@doc raw"""
    plot_surface_base!(x, y, Z; fig=nothing, ax=nothing, hm=nothing, xlabel=L"$x$", ylabel=L"$y$", zlabel=L"$z$", title="Title", xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, colormap=:jet)

Basic surface plot of `Z` vs `x` and `y`.
* `fig`: Figure
* `ax`: Axis of the figure
* `hm`: Colormap of the figure

You can overwrite on a plot by specifying `fig`, `ax` and `hm`. Otherwise, it will create a new plot. Return `fig`, `ax` and `hm`.

"""
function plot_surface_base!(x, y, Z; fig=nothing, ax=nothing, hm=nothing, xlabel=L"$x$", ylabel=L"$y$", zlabel=L"$z$", title="Title", xlims=nothing, ylims=nothing, zlims=nothing, xscale=log10, yscale=log10, zscale=log10, colormap=:jet)
  if fig === nothing || ax === nothing || hm === nothing
    if xlims === nothing
      xlims = (minimum(x), maximum(x))
    end

    if ylims === nothing
      ylims = (minimum(y), maximum(y))
    end

    if zlims === nothing
      zlims = (minimum(Z)/1.1, maximum(Z)*1.1)
    end

    #TODO: Workaround for log scale axes scales Axis3. Change this when Makie deals with it (https://docs.makie.org/stable/reference/blocks/axis3/).
    if xscale == log10
      xp = log10.(x)
      xticks = LinRange(log10(xlims[1]),log10(xlims[2]),4)
    else
      xp = x
      xticks = LinRange(xlims[1],xlims[2],4)
    end

    if yscale == log10
      yp = log10.(y)
      yticks = LinRange(log10(ylims[1]),log10(ylims[2]),4)
    else
      yp = y
      yticks = LinRange(ylims[1],ylims[2],4)
    end

    if zscale == log10
      Zp = log10.(Z)
      zticks = LinRange(log10(zlims[1]),log10(zlims[2]),4)
    else
      Zp = Z
      zticks = LinRange(zlims[1],zlims[2],4)
    end
    fig, ax, hm = surface(xp, yp, Zp; 
      axis=(;type=Axis3),
      colormap = colormap,
      colorrange = (zlims)
    )
    Colorbar(fig[1, 2], hm)
    resize!(fig.scene,(400,300))
    ax.xlabel=xlabel
    ax.xlabelsize=fontsizeAxisLabel
    ax.xticklabelsize=fontsizeAxis
    ax.ylabel=ylabel
    ax.ylabelsize=fontsizeAxisLabel
    ax.yticklabelsize=fontsizeAxis
    ax.zlabel=zlabel
    ax.zlabelsize=fontsizeAxisLabel
    ax.zticklabelsize=fontsizeAxis
    ax.title=title
    ax.titlesize=fontsizeTitle

    #TODO: Workaround for log scale axes scales Axis3. Change this when Makie deals with it (https://docs.makie.org/stable/reference/blocks/axis3/).
    if xscale == log10
      ax.xticks = (xticks, ["10" * Makie.UnicodeFun.to_superscript(round(Int64, v)) for v in xticks])
    end

    if zscale == log10
      ax.yticks = (yticks, ["10" * Makie.UnicodeFun.to_superscript(round(Int64, v)) for v in yticks])
    end
    
    if zscale == log10
      ax.zticks = (zticks, ["10" * Makie.UnicodeFun.to_superscript(round(Int64, v)) for v in zticks])
    end
    
    xlims!(ax,(xticks[1],xticks[end]))
    ylims!(ax,(yticks[1],yticks[end]))
    zlims!(ax,(zticks[1],zticks[end]))
  else
    empty!(ax.scene.plots)
    surface!(ax, x, y, Z; colormap=hm.colormap, colorrange=hm.colorrange)
  end
  return fig, ax, hm
end

@doc raw"""
    plot_2D_slices_base!(x, y, Z; dir=1, nslices=8, fig=nothing, ax=nothing, cb=nothing, slabel=L"$s$", zlabel=L"$z$", clabel=L"$c$", title="Title", slims=nothing, zlims=nothing, clims=nothing, sscale=log10, zscale=log10, cscale=log10, colormap=:inferno)  

Plot `nslices` slices of `Z` in direction `dir`.
* `fig`: Figure
* `ax`: Axis of the figure
* `cb`: Colormap of the figure

You can overwrite on a plot by specifying `fig`, `ax` and `cb`. Otherwise, it will create a new plot. Return `fig`, `ax` and `cb`.


"""
function plot_2D_slices_base!(x, y, Z; dir=1, nslices=8, fig=nothing, ax=nothing, cb=nothing, slabel=L"$s$", zlabel=L"$z$", clabel=L"$c$", title="Title", slims=nothing, zlims=nothing, clims=nothing, sscale=log10, zscale=log10, cscale=log10, colormap=:inferno)  
  if dir == 1
    s = y
    c = x
  elseif dir == 2
    s = x
    c = y
  else
    @error "dir must be 1 or 2. wave action with three arguments not yet implemented"
  end

  if fig === nothing || ax === nothing || cb === nothing
    if slims === nothing
      slims = (minimum(s), maximum(s))
    end

    if clims === nothing
      clims = (minimum(c), maximum(c))
    end

    if zlims === nothing
      if zscale == log10
        inds = Z .> 0
        zlims = (minimum(Z[inds])/1.1, maximum(Z[inds])*1.1)
      else
        zlims = (minimum(Z)/1.1, maximum(Z)*1.1)
      end
    end

    fig = Figure(size = (400, 300))  
    ax = Axis(fig[1, 1],
      xlabel = slabel,
      xlabelsize=fontsizeAxisLabel,
      xscale=sscale,
      xticklabelsize=fontsizeAxis,
      ylabel=zlabel,
      ylabelsize=fontsizeAxisLabel,
      yscale=zscale,
      yticklabelsize=fontsizeAxis,
      title=title,
      titlesize=fontsizeTitle,
    )
    xlims!(ax,slims)
    ylims!(ax,zlims)

    cb = Colorbar(fig[1, 2], limits = clims, colormap = colormap, label = clabel, scale = cscale)
  end

  empty!(ax)
  d = fld(length(s),nslices)
  for i in 1:d:length(c)
    if dir == 1
      z = Z[i,:]
    elseif dir ==2
      z = Z[:,i]
    else
      @error "dir must be 1 or 2. wave action with three arguments not yet implemented"
    end
    lines!(ax, s, z, colormap=cb.colormap, color=c[i], colorrange=cb.limits, colorscale=cb.scale)  
  end
  
  return fig, ax, cb
end

@doc raw"""
    plot_theo!(ax, kk::Vector{Float64}, nk::Vector{Float64}; color="red")

Add the theoretical line `nk` vs `kk` on figure with axis `ax`.

"""
function plot_theo!(ax, kk::Vector{Float64}, nk::Vector{Float64}; color="red")
  inds =  nk .> 0
  lines!(ax, kk[inds],nk[inds],linestyle=:dash,color=color)
end