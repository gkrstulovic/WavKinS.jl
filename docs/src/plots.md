# Plots

We define methods for plots that are wrappers of methods of the [`Makie.jl`](https://docs.makie.org/stable/) package. Many methods have optionnal arguments (e.g. `title`, `xlabel`, `xscale`, ...) consistent with [`Makie.jl`](https://docs.makie.org/stable/), so please refer to this well documented package for more informations. Wavkins plots allows one to plot standard quantities used in the theory (spectra of different quantities, associated fluxes, global quantities, etc.) quite simply. 

!!! warning 
    If ones tries to plot a spectrum with negative or zero values in a log-log plot, Makie crashes WavKinS. In the next versions, the plotting routines will be decoupled form the main code and this issue will be solved.


## Basics

```@autodocs
Modules = [WavKinS]
Order   = [:type, :function]
Pages   = ["src/plot/base.jl"]
```

## Spectra 

```@autodocs
Modules = [WavKinS]
Order   = [:type, :function]
Pages   = ["src/plot/spectra.jl"]
```

## Fluxes 

```@autodocs
Modules = [WavKinS]
Order   = [:type, :function]
Pages   = ["src/plot/fluxes.jl"]
```