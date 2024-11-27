# How to create a new solver

## The minimal steps to implement a new physical systems are

1. If none of the structures defined in `/src/WavKinS_structures.jl` can be used in your problem, create a new one. In this case, most of the routines should probably be adapted.
2. Create a new folder in `/src/physical_systems/`. It must contain:
  - `basics.jl`: definitions for your physical system (interaction coefficients, dispersion relation, parameterization of the resonant manifold, etc.).
  - `structure.jl`: fields and other variables needed for a typical simulation of your new system (diagnostics, working fields, wave-action, etc.).
  - `collision_integral.jl`: defines the collision integral `St_k!(Run)`.
  - `special.jl`: special diagnostics are defined here (optional).  
  If you think you're working on something generic (i.e. that can be used for other physical systems), maybe you should write it in different parts of the code (see code structure in `index.html`). If you're willing to work on a variation of an implemented physical system, you must write it in the same folder (e.g. `Acoustic2D` and `Acoustic3D` are both in `/src/physical_systems/Acoustic/`) for consistency.



## Some general advice for good practice

1. Keep consistent names of fields inside all structures. In doing so, existing routines can be used with your new simulation structures. For example, time-stepping routines in `/src/time_stepping/WavKinS_time_steping.jl` are generic.
2. Write test files to check what you have implemented for your physical system (parameterization of the resonant manifold, conservation of the dynamical invariant(s) by the collisional integral, ...) in the spirit of `/run/tests/physical_systems/tests_Acoustic2D.jl` and `/run/tests/physical_systems/tests_Petviashvilli.jl`.
3. If you have implemented more general tools (integration, interpolation, time stepping, ...), it should also be tested (see, e.g., `/run/tests/tests_integration_interpolation.jl` and `/run/tests/tests_misc.jl`) and documented.
4. Please give a short description of your physical system, how you compute the collision integral, and autodocumentation of methods and structures in `docs/src/physical_systems/`. You should add a link to this documentation in [`Currently available systems`](@ref "Currently available systems") subsection.