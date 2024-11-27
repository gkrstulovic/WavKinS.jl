# The structures and the constructors of the Stratified runs go here


@doc raw"""
    Stratified_Asymp

Simulation structure for internal gravity wave turbulence in the hydrostatic limit ``k_h \ll |k_z|``. It contains

    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments

    Nk::wave_spectrum_khkz #wave action
    Sk::wave_spectrum_khkz #collisional integral
    F1::wave_spectrum_khkz #working field
    kin_box::kinematic_box #kinematic box: meshes for a (p) and q
    FStp::Vector{Vector{Float64}} #Array of working fields for multithreading
    FStq::Vector{Vector{Float64}} #Array of working fields for multithreading
    partition::Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}} #partition for multithreading

    # Type of interpolation, integration, and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    integ_scheeme::Integ_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    N::Float64 #buoyancy (Brünt Väisälä) frequency
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation. 

"""
mutable struct Stratified_Asymp{Interp_Scheeme<:abstract_interpolation,Integ_Scheeme<:abstract_integration,Time_Stepping<:abstract_time_stepping}
    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments

    Nk::wave_spectrum_khkz #wave action
    Sk::wave_spectrum_khkz #collisional integral
    F1::wave_spectrum_khkz #working field
    kin_box::kinematic_box #kinematic box: meshes for a (p) and q
    FStp::Vector{Vector{Float64}} #Array of working fields for multithreading
    FStq::Vector{Vector{Float64}} #Array of working fields for multithreading
    partition::Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}} #partition for multithreading

    # Type of interpolation, integration, and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    integ_scheeme::Integ_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    N::Float64 #buoyancy (Brünt Väisälä) frequency
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation. 
end

@doc raw"""
    Stratified_Asymp(Nk::wave_spectrum_khkz; interp_scheeme=bilin_interp_khkz,integ_scheeme=integrate_with_log_bins_khkz,time_stepping_scheeme=RK2_step,Mp=0,Mq=0,nthh=0,nthz=0,N=1.0)

Constructor of a [`Stratified_Asymp`](@ref) structure. Optionally we set interpolation and time-stepping scheemes:

* `interp_scheeme`: [`bilin_interp_khkz`](@ref) (default), [`cpow_interp_khkz`](@ref)
* `integ_scheeme` : [`integrate_with_log_bins_khkz`](@ref) (default), [`integrate_with_cpow_khkz`](@ref)
* `time_stepping_scheeme` : [`Euler_step`](@ref), [`RK2_step`](@ref) (default), [`RK4_step`](@ref)

You can also change:

* `Mp`, `Mq`: number of grid points of the kinematic box in ``p`` and ``q`` directions
* `nthh`, `nthz`: number of slices in ``k_h`` and ``k_z`` for multithreading
* `N=1.0`: buoyancy frequency

"""
function Stratified_Asymp(Nk::wave_spectrum_khkz; interp_scheeme=bilin_interp_khkz,integ_scheeme=integrate_with_log_bins_khkz,time_stepping_scheeme=RK2_step,Mp=0,Mq=0,nthh=0,nthz=0,N=1.0)

    name = "Stratified_Asymp"
    Nk_arguments = 2
    ω = (kh, kz) -> ω_Stratified_Asymp(kh, kz; N)

    kkh = Nk.kkh
    khmin = kkh[1]
    khmax = kkh[end]
    Mh = Nk.Mh
    kkz = Nk.kkz
    kzmin = kkz[1]
    kzmax = kkz[end]
    Mz = Nk.Mz

    Stk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)

    # We define the kinematic box for computing the collision integral
    if Mq == 0
        Mq = 2*Mh
    end
    amin = (khmin/Mh)*ones(length(kkh))
    amax = kkh
    Ma = Int.(max.(8*ones(length(kkh)),1:Mh))
    qmin = khmin/Mq
    qmax = 2*khmax


    kin_box=kinematic_box(amin,amax,Ma,qmin,qmax,Mq)

    Mp = length(kin_box.aa[end])

    FStp = Vector{Vector{Float64}}(undef, Threads.nthreads())
    FStq = Vector{Vector{Float64}}(undef, Threads.nthreads())

    for ith = 1:Threads.nthreads()
        FStp[ith] = Vector{Float64}(undef,Mp)
        FStq[ith] = Vector{Float64}(undef,Mq)
    end

    F1 = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)

    diags = diagnostic_container(Mh,Mz)

    #Set partition for multithreading
    partition = get_partition(Mh,Mz,Threads.nthreads(),nthh,nthz)

    t = 0.0
    dimension = 2
    dΩ = 2 * pi

    FD = force_dissipation(Mh,Mz)

    interpScheeme = interp_scheeme(kkh,kkz)
    integScheeme = integ_scheeme()
    time_steppingScheeme = time_stepping_scheeme(Mh,Mz)

    Stratified_Asymp(name, Nk_arguments, ω, Nk, Stk, F1, kin_box, FStp, FStq, partition, interpScheeme, 
        integScheeme, time_steppingScheeme, diags, t, N, dimension, dΩ, FD)

end
