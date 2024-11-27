# Basics

We use structures for storing grids, dissipation coefficients, parameters, diagnostics, etc. For all physical systems, we define a general structure containing all variables and parameters necessary for a simulation (see e.g. [`Acoustic2D`](@ref "Acoustic2D")).

## Structures
```@autodocs
Modules = [WavKinS]
Order   = [:type, :function]
Pages   = ["src/WavKinS_structures.jl"]
```

## Parameters
```@autodocs
Modules = [WavKinS]
Order   = [:type, :function]
Pages   = ["src/WavKinS_parameters.jl"]
```

## Diagnostics
```@autodocs
Modules = [WavKinS]
Order   = [:type, :function]
Pages   = ["src/WavKinS_diagnostics.jl"]
```

## Outputs
```@autodocs
Modules = [WavKinS]
Order   = [:type, :function]
Pages   = ["src/WavKinS_outputs.jl"]
```