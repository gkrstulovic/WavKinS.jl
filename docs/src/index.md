# WavKinS: Wave Kinetic equations Solver

Wave kinetic equations arrive in many physical systems, such as plasmas, geophysical fluids, elastic plates, turbulent Bose-Einstein condensates, and many others. When the nonlinearity is small, the theory of Weak Wave Turbulence (WWT) provides a close equation describing the evolution of the wave spectrum, roughly describing the wave amplitude or energy at a given scale. WavKinS currently solves kinetic equations of the 3—and 4-wave types, as described below.

WavKinS is software developed in the [Julia](https://julialang.org) language. It aims to provide a simple and efficient general solver of wave kinetic equations. Thanks to the structures and routines already developed in WavKinS, we expect developing a new solver for a specific physical system to be simple.

!!! warning
    WavKinS is a project under development; many routines could still be improved for accuracy or efficiency. All remarks and contributions from the scientific community are welcome! 

## License

WavKinS.jl is an Open Source project distributed with [EUPL-1.2](https://interoperable-europe.ec.europa.eu/collection/eupl) license agreement. 

!!! info "Citing WavKinS and intellectual recognition"
    In addition to the EUPL obligations; the use of this Work for the analysis, creation or modification of data imposes on the Licensee, in each of their publications resulting from the exploitation of these data, the obligation to cite the first articles where this Work was used: [V. Labarre, G. Krstulovic and S. Nazarenko arXiv:2407.11469](https://arxiv.org/abs/2407.11469) and [Y. Zhu, G. Krstulovic and S. Nazarenko arXiv:2408.15163](https://arxiv.org/abs/2408.15163).


## General description

This program solves the Wave Kinetic Equation (WKE)

$$\frac{\mathrm{d}n_{\bf k}}{\mathrm{d}t} = St_{\bf k} + f_{\bf k} - d_{\bf k} n_{\bf k}$$  

where $n_{\bf k}$ is the wave-action spectrum, $St_{\bf k}$ is the collisional integral, $f_{\bf k}$ is the forcing, and $d_{\bf k}$ the dissipation coefficient. Their definitions depend on the physical problem. See e.g. [Nazarenko, Springer, Volume 825 of Lecture Notes in Physics (2011)](https://link.springer.com/book/10.1007/978-3-642-15942-8) for an introduction on the topic. 

The WKE can be classified by the number of waves interacting, which results in different degrees of nonlinearity and integrations to be performed. In the current version of the code, we consider two types of systems that can be generically written as follows.

### 3-waves canonical systems

In this case the collisional integral reads

```math
St_{\bf k}=\int (\mathcal{R}^{\bf k}_{12}-\mathcal{R}^{1}_{2 {\bf k}}-\mathcal{R}^{2}_{{\bf k}1})\mathrm{d}{\bf k}_1\mathrm{d}{\bf k}_2
```

with $\mathcal{R}^{\bf k}_{12}=|V_{{\bf k}12}|^2(n_1n_2-n_1n_{\bf k}-n_{\bf k}n_2)\delta({{\bf k}-{\bf k}_1-{\bf k}_2})\delta(\omega_{\bf k}- \omega_1- \omega_2)$, $V_{{\bf k}12}$ is the interaction coefficient, and $\omega_{\bf k}$ is the wave frequency. Systems with interactions involving more than 3 waves and/or having a noncanonical structure can be implemented as well.

Note that if wavevectors are in dimension $d$, then there are $2d - d - 1=d-1$ integrals to be performed for each value of a external wavevector ${\bf k}$. Therefore, if we use $M$ points to discretise each coordinate of the wavevector space, then the numerical cost will be generically of the order of $M^{2d -1}$ operations. However, assuming full isotropy or isotropy in a plane can drastically reduce the numerical cost, making calculations possible.

### 4-waves canonical systems

In these systems there is one extra integral
```math
St_{\bf k}=\int |T^{{\bf k 1}}_{\bf 2 3}|^2\delta({{\bf k}+{\bf k}_1-{\bf k}_2}-{\bf k}_3)\delta(\omega^{{\bf k 1}}_{\bf 2 3}) n_{\bf k}n_{\bf 1}n_{\bf 2}n_{\bf 3}\left( \frac{1}{n_{\bf k}} + \frac{1}{n_{\bf 1}} -\frac{1}{n_{\bf 2}} -\frac{1}{n_{\bf 3}} \right)\mathrm{d}{\bf k}_1\mathrm{d}{\bf k}_2\mathrm{d}{\bf k}_2
```
with $\omega^{{\bf k 1}}_{\bf 2 3}=\omega_{\bf k}+ \omega_{\bf 2}- \omega_{\bf 3}-\omega_{\bf 4}$ and $T^{{\bf k 1}}_{\bf 2 3}$ the interaction coefficient.

For these systems, there are $3d-d-1=2d-1$ integrals to be performed, which implies in general a numerical cost of the order of $M^{3d-1}$ operations.

### Methods

We mostly use logarithmic grids to span larger range of scales (wavevectors). The code includes standard interpolation and integration schemes adapted to these grids. Standard time evolution schemes (Euler, RK2, RK4, and Splitting methods) with adaptative time steps are also implemented. For most physical systems, we define basic outputs (e.g. wave-action, energy, and associated fluxes), that can be saved during the simulations.

Wavkins intends to be a generic and modular solver in which new physical systems should be easily implemented, see section [How to create a new solver](@ref "How to create a new solver"). 

!!! warning
    WavKinS is still under development. Although it has been carefully tested you can use at your own digression. If you find any bug or mathematical inconsistency please report it to make this solver


## Physical systems 

We have developed solvers for several physical systems. These solvers have been written in past and current scientific research projects for scientific purposes.

### Currently available systems

- [`Acoustic2D`](@ref "Acoustic2D solver"): Isotropic two-dimensional acoustic wave turbulence.

- [`Acoustic3D`](@ref "Acoustic3Dsolver"): Isotropic three-dimensional acoustic wave turbulence. 

- [`Bogoliubov3D`](@ref "Bogoliubov3Dsolver"): Isotropic three-dimensional Bogoliubov wave turbulence.

- [`MMT`](@ref "MMT solver"): Majda-McLaughlin-Tabak model.

- [`NLS3D`](@ref "NLS3Dsolver"): Isotropic three dimensional non-linear Schrödinger equation.

- [`Petviashvilli_Asymp`](@ref "Petviashvilli_Asympsolver"): Waves systems based on the Petviashvilli equation in the strongly anisotropic limit.

- [`Stratified_Asymp`](@ref "Stratified_Asympsolver"): Three-dimensional internal gravity wave turbulence in the strongly anisotropic (hydrostatic) limit. 

- [`Smoluchowski`](@ref "Smoluchowskisolver"): Smoluchowski coagulation equation.

### Under development

- `Acoustic2D_kxky`: two-dimensional acoustic wave turbulence with anisotropic forcing and dissipation.

- `IGW_Asymp` : internal and inertial gravity waves in the strongly anisotropic (hydrostatic) limit. 
- `IGW` :  non-hydrostatic internal and inertial gravity waves

## Code structure

The code is decomposed as follows:

- `src`: Source files of the code. It is divided into several parts:
    * `grid`: Structures and functions for the grids.
    * `integration`: Tools used for integrations, mainly for log grids. 
    * `interpolation`: Tools used for interpolations, mainly for log grids.
    * `misc`: Miscellaneous tools, that are still not classified.
    * `physical_systems`: Structures and necessary functions for the simulations of wave systems. We use a sub-folder for each systems, which contains: 
        - `basics.jl`: Definitions of the frequency, interaction coefficients, resonant manifold, ...
        - `collision_integral.jl`: Define the method [`St_k!`](@ref) for computing the collision integral. [`St_k!`](@ref) is overloaded, allowing to have the same syntax for all solvers. 
        - `structure.jl`: Structure for the run.
        - `special.jl`: Non-standard diagnostics/functions (not for all systems). 
        We give common definitions and structures, used for several wave systems, in [`basics_all.jl`]. For some wave systems, we packed several variations in the same sub-folders for consistency (e.g. `Acoustic` contains structures and definitions of isotropic and anisotropic Acoustic wave systems whose dispersion relation is the same, but not the spatial dimension and/or the forcing).
    * `plot`: Standard `WavKins` plots. They are all made with the [`Makie.jl`](https://docs.makie.org/stable/) package, using the `GLMakie.jl` backend. #TODO: Move this part elsewhere when separated from the main code.
    * `time_stepping`: Standard functions for time stepping.
    * `WavKinS_diagnostics.jl`: Definitions of standard diagnostics (total wave-action, total energy, total energy dissipation, energy spectrum, wave-action fluxes, energy fluxes, ...). 
    * `WavKinS_outputs.jl`: Inputs/Outputs functions. 
    * `WavKinS_parameters.jl`: Structures for simulation parameters (output directory, which outputs, final time, ...).
    * `WavKinS_partitions_threads.jl`: Tools for multithreading. In particular, the functions for getting partitions of nodes among threads.
    * `WavKinS_structures.jl`: Very standard structures used in WavKins simulations (spectra, params containers, ...).
    
- `run`: Run simulations and test the code.
    * `simple`: Simple simulations for various physical systems.  
    * `tests`: Test the different parts of the code (see the folder `src`).

- `postproc`: Example scripts (Julia, Python, and Matlab) that show how we can analyze outputs of some simulations.

- `materials`: Useful references and notebooks that explain the theory and some technical parts of the code.

- `docs`: Documentation.


## Contributors

Main:
* [Giorgio Krstulovic](https://gkrstulovic.gitlab.io/) (krstulovic@oca.eu) 
* [Vincent Labarre](https://scholar.google.com/citations?user=83UL2ZsAAAAJ&hl=en) (vincent.labarre@polytechnique.edu) 

Other developers:
* [Ying Zhu](https://scholar.google.fr/citations?user=N43KuycAAAAJ&hl=en)
* [Guillaume Costa](https://scholar.google.com/citations?user=iSn_fuYAAAAJ&hl=fr)

Scientific collaborators:
* [Sergey Nazarenko](https://scholar.google.com/citations?user=EPW6UlQAAAAJ&hl=en)
* [Juan Ignacio Polanco](https://jipolanco.gitlab.io/)
