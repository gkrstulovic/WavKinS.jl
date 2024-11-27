# Computing a collisional integral

> Look at `README.md`, and install the code first.

In this simple tutorial we show how to initialise the WavKinS, to define a couple of basics structures needed to make calculations, set an initial waveaction and finally compute the collisional integral (the right hand side of the WKE).

!!! note
    Before running any of the commands below you need to indicate to Julia the path to the folder where WavKinS has been installed, here denoted as `PATH_TO_WAVKINS_FOLDER`. In Linux bash, this is set by setting the `JULIA_LOAD_PATH` environment variable

    ```bash
    export JULIA_LOAD_PATH="PATH_TO_WAVKINS_FOLDER:$JULIA_LOAD_PATH" 
    ```

## Defining the basic structures.

One of the most basic structures in WavKinS is the [`wave_spectrum`](@ref) for isotropic systems. The constructor requires to set the number of nodes, and the smallest and the largest wavenumbers. To initialise a [`wave_spectrum`](@ref) we do the following

```julia
## Defining the wave_spectrum structure
using WavKinS
M = 1024
kmin = 1e-3
kmax = 1e0
Nk = wave_spectrum(kmin, kmax, M)
``` 

We have now the waveaction `Nk.nk` initialised with zeros and the logarithmic grid `Nk.kk`. We are going set the waveaction to the Rayleigh–Jeans distribution
```math
n_{\bf k} = \frac{1}{\omega_{\bf k}}
```
This distribution is a trivial solution of the WKE, which means that $St_{\bf k}=0$ ([Nazarenko, Springer, Volume 825 of Lecture Notes in Physics (2011)](https://link.springer.com/book/10.1007/978-3-642-15942-8)). Here $\omega_{\bf k}$ is the wave dispersion relation that depends on each physical system.

For this tutorial, we will use the [`Acoustic2D`](@ref) solver. For each solver, there is one structure containing all the methods, arrays, and parameters needed for a simulations. We initialise it as follows

```julia
Run = Acoustic2D(Nk)
```
Here we did not provide any optional parameters. If desired, one can set the time-stepping, interpolation and integration method, physical parameters etc. We can now set our Rayleigh–Jeans distribution using the dispersion relation included in `Run`
```julia
kk = Run.Nk.kk # we rename the wavevector mesh for shorter use
nk = Run.Nk.nk # we rename the waveaction mesh for shorter use
@. nk = 1.  / Run.ω(kk) # the @. is the julia macro for pointwise operations
```
Finally, we compute the collisional integral
```julia
St_k!(Run) # compute the collision integral
```
The collisional integral is stored in the [`wave_spectrum`](@ref)  structure `Run.Sk`, and therefore, its values are accessible in `Run.Sk.nk`.


To check the error of the collisional integral, we compute the $L^2$ norm $||St ||=\sqrt{\int_0^\infty 2\pi k |St_k|^2 dk}$

```julia
 @. Run.Sk.nk = 2pi * kk * Run.Sk.nk^2 # we use the same waveaction structure to compute the square

L2norm = sqrt(integrate_with_log_bins(Run.Sk)); # we integrate collisional integral squared using WavKinS integration
println("The L2-norm is ",L2norm)
```
which produces the output
```julia-repl
The L2-norm is 0.0005670508265665057
```

!!! info
    Note that the Rayleigh–Jeans distributions are not integrable functions, however all integrals in WavKinS have sharp cut-off at `kmax`, giving mathematical sense to these solutions.

## Running this example

To run the above example for increasing values of the number of points `M`, just copy the following lines to a file (here called `computing_collisional_integral_tutorial.jl`)
```julia
using WavKinS

## Defining the wave_spectrum structure

kmin = 1e-3;
kmax = 1e0;

for M ∈ 2 .^ (5:12)
    Nk = wave_spectrum(kmin, kmax, M)

    # Creating the Run Acoustic2D, filling the waveaction and computing the collisional integral
    Run = Acoustic2D(Nk)
    kk = Run.Nk.kk # we rename the wavevector mesh for shorter use
    nk = Run.Nk.nk # we rename the wavefunction mesh for shorter use
    @. nk = 1.0 / Run.ω(kk) # the @. is the julia macro for pointwise operations


    # Compute the collisional integral and print th
    St_k!(Run) # compute the collision integral

    # Compute the L2-norm of the collisional integral
    @. Run.Sk.nk = 2pi * kk * Run.Sk.nk^2 # we use the same waveaction structure to compute the square
    L2norm = sqrt(integrate_with_log_bins(Run.Sk)) # we integrate collisional integral squared using WavKinS integration

    println("M =", M, " L2-norm = ", L2norm)

end
```

Then to run the file with 4 threads, simply execute from the terminal

```bash
julia --threads 4 computing_collisional_integral_tutorial.jl
```

which produces
```julia_repl
M =32 L2-norm = 0.11497899775058602
M =64 L2-norm = 0.039508648493920724
M =128 L2-norm = 0.013842127550835866
M =256 L2-norm = 0.00486925111983586
M =512 L2-norm = 0.0016955053736738902
M =1024 L2-norm = 0.0005670508265665057
M =2048 L2-norm = 0.00016406073037112257
M =4096 L2-norm = 2.4680754912948397e-5
```

We observe that error decreases as $M^{-2}$.