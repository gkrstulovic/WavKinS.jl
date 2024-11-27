# Defining a new diagnostic

Suppose that you would like to run WavKinS and compute a new diagnostic that has not been included in the original physical system. You would like also that WavKinS writes to the disk this new measurement together the default ones. 

In this tutorial, we assume that you know how to define the basic WavKinS structures. We now explain how to add a new diagnostic computed from the waveaction $n_{\bf k}$.


!!! note
    Before running any of the commands below you need to indicate to Julia the path to the folder where WavKinS has been installed, here denoted as `PATH_TO_WAVKINS_FOLDER`. In Linux bash, this is set by setting the `JULIA_LOAD_PATH` environment variable

    ```bash
    export JULIA_LOAD_PATH="PATH_TO_WAVKINS_FOLDER:$JULIA_LOAD_PATH" 
    ```

## Basic structures for diagnostic

For simplicity, in this tutorial we also consider the [`Acoustic2D`](@ref) physical system. As in the previous tutorial, we initialise the waveaction and [`Acoustic2D`](@ref) structures and fill the wave action with something not trivial (the Rayleigh–Jeans distribution). 

```julia
using WavKinS
M = 1024
kmin = 1e-3
kmax = 1e0
Nk = wave_spectrum(kmin, kmax, M)

Run = Acoustic2D(Nk)

kk = Run.Nk.kk # we rename the wavevector mesh for shorter use
nk = Run.Nk.nk # we rename the waveaction mesh for shorter use
@. nk = 1.  / Run.ω(kk) # the @. is the julia macro for pointwise operations
```
The default diagnostics are stored in `Run`, and initially empty or set to zeroes. They are stored as a dictionary of [`WavKinS.global_ouput`](@ref) and [`WavKinS.spectral_output`](@ref) (follow links for a precise definition):

```julia-repl
julia> Run.diags.glob_diag
Dict{String, WavKinS.global_ouput} with 4 entries:
  "N"     => global_ouput("Total wave action", Float64[])
  "Times" => global_ouput("Times of global quantities", Float64[])
  "H"     => global_ouput("Total energy", Float64[])
  "Disp"  => global_ouput("Total dissipation", Float64[])
```
and
```julia-repl
julia> Run.diags.sp_outs
Dict{String, WavKinS.spectral_output} with 3 entries:
  "nk" => spectral_output("Wave action spectrum", [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], true)
  "Ek" => spectral_output("Energy spectrum", [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], true)
  "Pk" => spectral_output("Energy flux", [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0  …  0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], true)
```

We want to add two new diagnostics, for instance a global quantity $B$ and a spectral one $B_k$
```math
B= \int \rho_k n_k d{\bf k}, \quad   B_k=\rho_k n_k
```
where $\rho_k$ is some density function of $k$ only. Note that the invariant $B$ is computed over all the Fourier space. For the [`Acoustic2D`](@ref), $d{\bf k}$ can be replaced by $2\pi d{k}$.

Let us first to create a place where to store these two diagnostics using the constructors of the diagnostic and the flexibility of the dictionary

```julia
Run.diags.glob_diag["B"]= WavKinS.global_ouput("My new global diagnostic")
Run.diags.sp_outs["Bk"] = WavKinS.spectral_output("My new spectral diagnostic",Nk.M)
```
For instance, [`WavKinS.global_ouput`](@ref) contains now
```julia-repl
julia> Run.diags.glob_diag
Dict{String, WavKinS.global_ouput} with 5 entries:
  "B"     => global_ouput("My new global diagnostic", Float64[])
  "N"     => global_ouput("Total wave action", Float64[])
  "Times" => global_ouput("Times of global quantities", Float64[])
  "H"     => global_ouput("Total energy", Float64[])
  "Disp"  => global_ouput("Total dissipation", Float64[])
```
We see now that $B$ is included in the global diagnostics.

## Computing the new diagnostic

Let's first define the density $\rho_k$ that we will use for the diagnostics (an arbitrary function of $k$)

```julia
function ρk(k)
    return k^4 /(1 + k^2)
end
```
We are now going to compute the default and the new diagnostics. First we call the WavKinS routines for the default diagnostics

```julia
get_global_diagnostics!(Run)
compute_spectral!(Run)
```

Finally, we compute the new diagnostics and add them to diagnostics
```julia
# For the global diagnostic
B = WavKinS.total_integral_density(Run, ρk) # we use WavKinS routine to compute B
push!(Run.diags.glob_diag["B"].out, B) # we append it to the diagnostic

# For the spectral diagnostic
Bk = Run.diags.sp_outs["Bk"].sp;
nk = Run.Nk.nk;
kk = Run.Nk.kk;
@. Bk = ρk(kk) * nk;  
```

The diagnostics have now
```julia-repl
julia> Run.diags.glob_diag
Dict{String, WavKinS.global_ouput} with 5 entries:
  "B"     => global_ouput("My new global diagnostic", [0.74606])
  "N"     => global_ouput("Total wave action", [6.28321])
  "Times" => global_ouput("Times of global quantities", [0.0])
  "H"     => global_ouput("Total energy", [3.14164])
  "Disp"  => global_ouput("Total dissipation", [0.0])
```
and
```julia-repl
julia> Run.diags.sp_outs
Dict{String, WavKinS.spectral_output} with 4 entries:
  "nk" => spectral_output("Wave action spectrum", [1000.0, 993.27, 986.586, 979.946, 973.352, 966.801, 960.295, 953.833, 947.414, 941.038  …  1.06266, 1.05551, 1.0484, 1.04135, 1.03434, 1.02738, 1.02046, 1.0136, 1.00678, 1.0], true)
  "Ek" => spectral_output("Energy spectrum", [0.00628319, 0.00632576, 0.00636861, 0.00641176, 0.00645521, 0.00649894, 0.00654297, 0.0065873, 0.00663194, 0.00667687  …  5.91271, 5.95277, 5.99311, 6.03371, 6.07459, 6.11575, 6.15719, 6.1989, 6.24…
  "Pk" => spectral_output("Energy flux", [3.47654e-11, 3.5239e-11, 3.57211e-11, 3.62117e-11, 3.6711e-11, 3.7219e-11, 3.77361e-11, 3.82622e-11, 3.87975e-11, 3.93423e-11  …  -9.75031e-7, -1.16121e-6, -1.35282e-6, -1.54999e-6, -1.75284e-6, -1.961…
  "Bk" => spectral_output("My new spectral diagnostic", [9.99999e-10, 1.02046e-9, 1.04135e-9, 1.06266e-9, 1.0844e-9, 1.10659e-9, 1.12924e-9, 1.15235e-9, 1.17593e-9, 1.19999e-9  …  0.44196, 0.448142, 0.454391, 0.460705, 0.467086, 0.473534, 0.48…
```

## Writing the diagnostic into a file

WavKinS uses the NCDF format using the [NCDataset](https://alexander-barth.github.io/NCDatasets.jl/stable/) Julia Package to mange the I/O. Once the new diagnostics are created and stored in [`WavKinS.global_ouput`](@ref) and [`WavKinS.spectral_output`](@ref) they will be written automatically whenever the routines [`output_spectra`](@ref) and [`output_global`](@ref) are called.

To create the output we do the following

```julia
outputDir = "./"
init_IO(Run, outputDir) # Create an empty NCDF file
output_global(Run, outputDir) 
output_spectra(Run, outputDir)
```
After calling these routines, a file called in this case `WKE_Acoustic2D_data.nc` will created and written in `outputDir`. The name of the file depends on the solver used. To list the content of the file, we can simple make


```julia-repl
julia> using NCDataset
julia> file=NCDataset(outputDir * "WKE_" * Run.name * "_data.nc")
Dataset: WKE_Acoustic2D_data.nc
Group: /

Dimensions
   k = 1024
   timesSP = 1
   times = 1

Global attributes
  Dataset              = Acoustic2D
  Created with git commit id = unknown
Groups
  Dataset: WKE_Acoustic2D_data.nc
  Group: Global

  Variables
    B     (1)
      Datatype:    Float64 (Float64)
      Dimensions:  times
      Attributes:
       long_name            = My new global diagnostic

    N     (1)
      Datatype:    Float64 (Float64)
      Dimensions:  times
      Attributes:
       long_name            = Total wave action

    Times     (1)
      Datatype:    Float64 (Float64)
      Dimensions:  times
      Attributes:
       long_name            = Times of global quantities

    H     (1)
      Datatype:    Float64 (Float64)
      Dimensions:  times
      Attributes:
       long_name            = Total energy

    Disp     (1)
      Datatype:    Float64 (Float64)
      Dimensions:  times
      Attributes:
       long_name            = Total dissipation

  Dataset: WKE_Acoustic2D_data.nc
  Group: Spectral

  Variables
    k     (1024)
      Datatype:    Float64 (Float64)
      Dimensions:  k
      Attributes:
       long_name            = wave vectors

    TimesSP     (1)
      Datatype:    Float64 (Float64)
      Dimensions:  timesSP
      Attributes:
       long_name            = Times of spectral outputs

    nk     (1024 × 1)
      Datatype:    Float64 (Float64)
      Dimensions:  k × timesSP
      Attributes:
       long_name            = Wave action spectrum

    Ek     (1024 × 1)
      Datatype:    Float64 (Float64)
      Dimensions:  k × timesSP
      Attributes:
       long_name            = Energy spectrum

    Pk     (1024 × 1)
      Datatype:    Float64 (Float64)
      Dimensions:  k × timesSP
      Attributes:
       long_name            = Energy flux

    Bk     (1024 × 1)
      Datatype:    Float64 (Float64)
      Dimensions:  k × timesSP
      Attributes:
       long_name            = My new spectral diagnostic



julia> close(file)
closed Dataset
```
