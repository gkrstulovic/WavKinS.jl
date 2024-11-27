# The structures and the constructors of the Smoluchowski runs go here

@doc raw"""
    Smoluchowski

Simulation structure for Smoluchowski. It contains

    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Particle size
    K # Kernel

    Nk::wave_spectrum #wave action
    Sk::wave_spectrum #collisional integral
    F1::wave_spectrum #working field
    FSt::Vector{wave_spectrum} #Array of working fields for multithreading
    partition::Vector{UnitRange{Int64}} #partition for multithreading

    # Type of interpolation and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation. 

"""
mutable struct Smoluchowski{Interp_Scheeme<:abstract_interpolation,Time_Stepping<:abstract_time_stepping}
    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Particle size
    K # Kernel

    Nk::wave_spectrum #wave action
    Sk::wave_spectrum #collisional integral
    F1::wave_spectrum #working field
    FSt::Vector{wave_spectrum} #Array of working fields for multithreading
    partition::Vector{UnitRange{Int64}} #partition for multithreading

    # Type of interpolation and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation. 
end

@doc raw"""
    Smoluchowski(Nk::wave_spectrum; K=K_Smoluchowski_one, interp_scheeme=lin_interp,time_stepping_scheeme=RK2_step

Constructor of a [`Smoluchowski`](@ref) structure. Optionally we set interpolation and time-stepping scheemes:

* `interp_scheeme`: [`lin_interp`](@ref) (default), [`powexp_interp`](@ref), [`powGauss_interp`](@ref), [`BS_interp`](@ref)
* `time_stepping_scheeme`: [`Euler_step`](@ref), [`RK2_step`](@ref) (default), [`RK4_step`](@ref)

The optional parameter `K` is the coagulation kernel.

"""
function Smoluchowski(Nk::wave_spectrum; K=K_Smoluchowski_one, interp_scheeme=lin_interp,time_stepping_scheeme=RK2_step)

    name = "Smoluchowski"
    Nk_arguments = 1

    kk = Nk.kk
    kmin = kk[1]
    kmax = kk[end]
    M = Nk.M

    Stk = wave_spectrum(kmin,kmax,M)
    F1 = wave_spectrum(kmin,kmax,M)
    FSt = Vector{wave_spectrum}(undef, Threads.nthreads())
    for ith = 1:Threads.nthreads()
        FSt[ith] = wave_spectrum(kmin,kmax,M)
    end

    diags = diagnostic_container(M)

    #Set partition for multithreading
    partition = get_partition(M, Threads.nthreads())

    t = 0.0
    dimension = 1
    dΩ = 1.0

    FD = force_dissipation(M)

    interpScheeme = interp_scheeme(kk)
    time_steppingScheeme  = time_stepping_scheeme(M)

    Smoluchowski(name,Nk_arguments,ω_Smoluchowski,K,Nk, Stk, F1, FSt, partition,interpScheeme,time_steppingScheeme,diags,t,dimension,dΩ,FD)
end
