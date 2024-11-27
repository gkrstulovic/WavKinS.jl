[![documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://gkrstulovic.github.io/WavKinS.jl/dev/)

# WavKinS

Wave kinetic equations arrive in many physical systems, such as plasmas, geophysical fluids, elastic plates, turbulent Bose-Einstein condensates, and many others. When the nonlinearity is small, the theory of Weak Wave Turbulence (WWT) provides a close equation describing the evolution of the wave spectrum, roughly describing the wave amplitude or energy at a given scale. WavKinS currently solves kinetic equations of the 3—and 4-wave types, as described below.

WavKinS is software developed in the [Julia](https://julialang.org) language. It aims to provide a simple and efficient general solver of wave kinetic equations. Thanks to the structures and routines already developed in WavKinS, we expect developing a new solver for a specific physical system to be simple.

## License

WavKinS is an Open Source project distributed with [EUPL-1.2](https://interoperable-europe.ec.europa.eu/collection/eupl) license agreement. 

**Citing WavKinS and intellectual recognition:** In addition to the EUPL obligations; the use of this Work for the analysis, creation or modification of data imposes on the Licensee, in each of their publications resulting from the exploitation of these data, the obligation to cite the first articles where this Work was used: [V. Labarre, G. Krstulovic and S. Nazarenko arXiv:2407.11469](https://arxiv.org/abs/2407.11469) and [Y. Zhu, G. Krstulovic and S. Nazarenko arXiv:2408.15163](https://arxiv.org/abs/2408.15163).


## Installation

Download the Julia Programming Language https://julialang.org and follow the instructions for installation.

To install WavKinS open `julia` in project mode. From the `WavKinS` folder type:

```
> julia --project
```

then enter to the package mode by typing `]` and do
```
> instantiate
```

`Julia` should download and install all the necessary packages.

If you use `Visual Studio Code` for coding, it is useful to install the `julia` package [Revise.jl](https://timholy.github.io/Revise.jl/stable/).


## Generating Documentation

To generate a local copy of the documentation, to `/docs` folder and run

```
julia make.jl
```

The package [Documenter.jl](https://juliadocs.github.io/Documenter.jl/stable/) will generate a webpage in the folder `/docs/build`. After building, open the file `/docs/build/index.html`



## Getting started and running the code

The next tutorials show with examples the basic use of WavKinS. In all cases, we need to import the module `WavKinS`. Assuming that Julia knows where WavKinS has been installed, it is imported by doing
```julia
using WavKinS
```
before any call to WavKinS routines.


Before lunching Julia, you need to indicate the path to the folder where WavKinS has been installed, here denoted as `PATH_TO_WAVKINS_FOLDER`. In Linux bash, this is set by setting the `JULIA_LOAD_PATH` environment variable

```bash
export JULIA_LOAD_PATH="PATH_TO_WAVKINS_FOLDER:$JULIA_LOAD_PATH" 
```
