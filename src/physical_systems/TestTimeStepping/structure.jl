# The structures and the constructors of the Test system go here

@doc raw"""
    TestTimeStepping

Simulation structure for TestTimeStepping system. It contains

    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    
    Nk::wave_spectrum # wave action
    Sk::wave_spectrum # collisional integral
    F1::wave_spectrum # working field
    partition::Vector{UnitRange{Int64}} # partition for multithreading

    time_stepping::Time_Stepping

    t::Float64 # current time
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation. 

"""
mutable struct TestTimeStepping{Time_Stepping<:abstract_time_stepping}
    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    
    Nk::wave_spectrum # wave action
    Sk::wave_spectrum # collisional integral
    F1::wave_spectrum # working field
    partition::Vector{UnitRange{Int64}} # partition for multithreading

    time_stepping::Time_Stepping

    t::Float64 # current time
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation.  
end

@doc raw"""
    TestTimeStepping(Nk::wave_spectrum; time_stepping_scheeme=RK2_step)

Constructor of a [`TestTimeStepping`](@ref "TestTimeStepping") structure. 

* `Nk`: wave spectrum structure
* `time_stepping_scheeme`: [`Euler_step`](@ref "Euler_step"), [`RK2_step`](@ref "RK2_step") (default), [`ETD2_step`](@ref "ETD2_step"), ...

"""
function TestTimeStepping(Nk::wave_spectrum; time_stepping_scheeme=RK2_step)

    name = "TestTimeStepping"
    
    kk = Nk.kk
    kmin = kk[1]
    kmax = kk[end]
    M = Nk.M

    Stk = wave_spectrum(kmin,kmax,M)
    F1 = wave_spectrum(kmin,kmax,M)

    partition = get_partition(M, Threads.nthreads())

    t = 0.0

    FD = force_dissipation(M)

    time_steppingScheeme  = time_stepping_scheeme(M)

    TestTimeStepping(name,1,Nk,Stk,F1,partition,time_steppingScheeme,t,1,1.0,FD)
end
