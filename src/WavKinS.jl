module WavKinS
#   "WavKinS_structures.jl"
export field_grid_1D, wave_spectrum, wave_spectrum_khkz, kinematic_box

#   "WavKinS_basics.jl"
export update_coeff_interp!
export val_nk
export integrate, integrate_with_log_bins, integrate_quad_with_log_bins, integrate_with_grid, integrate_with_log_bins_khkz, integrate_with_cpow_khkz



export TestTimeStepping, Acoustic2D, Acoustic3D, Acoustic2D_khkz, Petviashvilli, Petviashvilli_Asymp, MMT, Stratified_Asymp, Bogoliubov3D, NLS3D, Smoluchowski



export simulation_parameters

#   "collision_integrals.jl"
export St_k!

#   "physical_systems.jl"

export total_waveaction

#   "drive_wavekin.jl"
#export eulerstep!
export Euler_step, RK2_step, RK4_step 
export advance!, get_T_nonlinear, adaptative_time_step
export init_temporal_scheeme!

#export ETD2coeff
#export ETDRK2step!

#   "WavKinS_diagnostics.jl"
export get_global_diagnostics!,store_spectral!,compute_spectral!
export energy, energy_flux!, waveaction_flux!, energy_dissipation, energy_injection, density_spectrum!, energy_spectrum!, isotropic_density_spectrum!, kh_density_spectrum!, kz_density_spectrum!, density_flux!, density_flux_isotropic!, density_flux_angular!, total_density_flux, total_density_abs_flux, energy_conservation_ratio, one

export output_global, output_spectra, init_IO, load_spectrum, load_spectral_for_change_mesh!

export change_mesh!, meshgrid, LogRange, area_ratio, area_ratio_grid, area_ratio_logbins, NewtonRaphson, NewtonRaphson_with_derivative

include("WavKinS_structures.jl")
include("WavKinS_diagnostics.jl")
include("WavKinS_outputs.jl")
include("WavKinS_parameters.jl")
include("WavKinS_partitions_threads.jl")


include("grid/WavKinS_grid.jl")

include("integration/WavKinS_integration_structs.jl")
include("integration/WavKinS_integration.jl")

include("interpolation/WavKinS_interpolation_structs.jl")
include("interpolation/WavKinS_interpolation.jl")

include("time_stepping/WavKinS_time_stepping_structs.jl")
include("time_stepping/WavKinS_time_stepping.jl")


# Plots
include("plot/base.jl")
include("plot/spectra.jl")
include("plot/fluxes.jl")
export plot_theo!, plot_1D_base!, plot_heatmap_base!, plot_surface_base!, plot_quiver_base, plot_2D_flux
export plot_density!, plot_wave_action!, plot_energy!, plot_1D_spectra_density!, plot_density_flux!, plot_wave_action_flux!, plot_energy_flux!, plot_2D_slices_base!, plot_2D_slices_khkz!, plot_2D_slices, plot_density_flux_isotropic!, plot_wave_action_flux_isotropic!, plot_energy_flux_isotropic!, plot_density_flux_angular!, plot_wave_action_flux_angular!, plot_energy_flux_angular!


include("misc/WavKinS_misc.jl")

include("physical_systems/basics_all.jl")


# Module TestTimeStepping
include("physical_systems/TestTimeStepping/basics.jl")
include("physical_systems/TestTimeStepping/structure.jl")
include("physical_systems/TestTimeStepping/collision_integral.jl")

# Module Acoustic 
include("physical_systems/Acoustic/basics.jl")
include("physical_systems/Acoustic/structure.jl")
include("physical_systems/Acoustic/collision_integral.jl")

# Module Bogoliubov 
include("physical_systems/Bogoliubov/basics.jl")
include("physical_systems/Bogoliubov/structure.jl")
include("physical_systems/Bogoliubov/collision_integral.jl")


# Module Petviashvilli
include("physical_systems/Petviashvilli/basics.jl")
include("physical_systems/Petviashvilli/structure.jl")
include("physical_systems/Petviashvilli/collision_integral.jl")
include("physical_systems/Petviashvilli/diagnostics.jl")

# Module MMT
include("physical_systems/MMT/basics.jl")
include("physical_systems/MMT/structure.jl")
include("physical_systems/MMT/collision_integral.jl")
include("physical_systems/MMT/diagnostics.jl")

# Module Stratified
include("physical_systems/Stratified/basics.jl")
include("physical_systems/Stratified/structure.jl")
include("physical_systems/Stratified/collision_integral.jl")
include("physical_systems/Stratified/diagnostics.jl")
include("physical_systems/Stratified/special.jl")
export compute_ωkz_spectrum, compute_ξkz_spectrum, plot_ωkz_spectrum!, plot_ξkz_spectrum!, plot_slices_ωkz_spectrum!, plot_slices_kzω_spectrum!, plot_slices_ξkz_spectrum!, plot_slices_kzξ_spectrum!, plot_ω_spectrum_density!

# Module NLS3D
include("physical_systems/NLS3D/basics.jl")
include("physical_systems/NLS3D/structure.jl")
include("physical_systems/NLS3D/collision_integral.jl")
include("physical_systems/NLS3D/diagnostics.jl")


# Module Smoluchowski
include("physical_systems/Smoluchowski/basics.jl")
include("physical_systems/Smoluchowski/structure.jl")
include("physical_systems/Smoluchowski/collision_integral.jl")

function __init__()
    PathToJulia=functionloc(WavKinS.eval)
    PathToJulia= PathToJulia[1]
    PathToJulia=PathToJulia[1:(end-10)]
    TTT=`git -C $PathToJulia rev-parse --short HEAD`
    commit = read(TTT,String)
    return commit[1:(end-1)]
end
commit_version = __init__()

end # module WavKinS
