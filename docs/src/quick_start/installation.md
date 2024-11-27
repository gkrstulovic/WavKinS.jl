# Installing WavKinS

The following steps are necessary to install WavKinS.

## Installation of Julia

Download the Julia Programming Language [https://julialang.org](https://julialang.org) and follow the instructions for installation.


## Installing WavkinS

To install WavKinS open `Julia` in project mode. From the `WavKinS` folder type:

```bash
> julia --project
```

then enter to the package mode by typing `]` and do
```julia-repl
(@v1.10) pkg> instantiate
```

`Julia` should download and install all the necessary packages.

If you use `Visual Studio Code` for coding, it is useful to install the `julia` package [Revise.jl](https://timholy.github.io/Revise.jl/stable/).


## Generating Documentation

To generate a local copy fo the documentation, go to  `/docs` folder and run

```bash
julia --project make.jl
```

The package [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/) will generate a webpage in the folder `/docs/build`. After building, open the file `/docs/build/index.html`



## Getting started and running the code

The next tutorials show with examples the basic use of WavKinS. In all cases, we need to import the module `WavKinS`. Assuming that Julia knows where WavKinS has been installed, it is imported by doing
```julia
using WavKinS
```
before any call to WavKinS routines.



!!! note
    Before lunching Julia, you need to indicate the path to the folder where WavKinS has been installed, here denoted as `PATH_TO_WAVKINS_FOLDER`. In Linux bash, this is set by setting the `JULIA_LOAD_PATH` environment variable

    ```bash
    export JULIA_LOAD_PATH="PATH_TO_WAVKINS_FOLDER:$JULIA_LOAD_PATH" 
    ```

