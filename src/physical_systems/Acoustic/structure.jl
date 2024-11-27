# The structures and the constructors of the Acoustic runs go here

@doc raw"""
    Acoustic2D

Simulation structure for Acoustic 2D wave turbulence. It contains

    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments

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
    a::Float64 # dispersive length
    c::Float64 # speed of sound
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation. 

"""
mutable struct Acoustic2D{Interp_Scheeme<:abstract_interpolation,Time_Stepping<:abstract_time_stepping}
    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments

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
    a::Float64 # dispersive length
    c::Float64 # speed of sound
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation 
end

@doc raw"""
    Acoustic2D(Nk::wave_spectrum; a=0.5, c=1.0,interp_scheeme=lin_interp,time_stepping_scheeme=RK2_step)

Constructor of a [`Acoustic2D`](@ref) structure. Optionally we set interpolation and time-stepping scheemes:

* `interp_scheeme`: [`lin_interp`](@ref) (default), [`powexp_interp`](@ref), [`powGauss_interp`](@ref), [`BS_interp`](@ref)
* `time_stepping_scheeme`: [`Euler_step`](@ref), [`RK2_step`](@ref) (default), [`RK4_step`](@ref)

The optional parameter `a=0.5` and `c=1.0` are the dispersive length and the speed of sound

"""
function Acoustic2D(Nk::wave_spectrum; a=0.5, c=1.0,interp_scheeme=lin_interp,time_stepping_scheeme=RK2_step)

    name = "Acoustic2D"
    Nk_arguments = 1
    ω = ω_Acoustic

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
    dimension = 2
    dΩ = 2*pi

    FD = force_dissipation(M)

    interpScheeme = interp_scheeme(kk)
    time_steppingScheeme  = time_stepping_scheeme(M)

    Acoustic2D(name,Nk_arguments,ω,Nk, Stk, F1, FSt, partition,interpScheeme,time_steppingScheeme,diags,t,a,c,dimension,dΩ,FD)
end

@doc raw"""
    Acoustic3D

Simulation structure for Acoustic 3D wave turbulence. It contains

    name::string #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments

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
    c::Float64 # speed of sound
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation. 


"""
mutable struct Acoustic3D{Interp_Scheeme<:abstract_interpolation,Time_Stepping<:abstract_time_stepping}
    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``. 1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments
   
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
    c::Float64 # speed of sound
    dimension::Int # physical dimension of the system (or of the isotropic sector)
    dΩ::Float64 # surface of the unit sphere (or of the isotropic sector)

    FD::force_dissipation # Contains all the terms about force and dissipation 
end

@doc raw"""
    Acoustic3D(Nk::wave_spectrum; a=0.5, c=1.0,interp_scheeme=lin_interp,time_stepping_scheeme=RK2_step

Constructor of a [`Acoustic3D`](@ref) structure. Optionally we set interpolation and time-stepping scheemes:

* `interp_scheeme`: [`lin_interp`](@ref) (default), [`powexp_interp`](@ref), [`powGauss_interp`](@ref), [`BS_interp`](@ref)
* `time_stepping_scheeme`: [`Euler_step`](@ref), [`RK2_step`](@ref) (default), [`RK4_step`](@ref)

The parameter `c=1.0` is the speed of sound

"""
function Acoustic3D(Nk::wave_spectrum; c=1.0,interp_scheeme=lin_interp,time_stepping_scheeme=RK2_step)

    name = "Acoustic3D"
    Nk_arguments = 1
    ω = ω_Acoustic

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

    diags = diagnostic_container(Nk.M)

    #Set partition for multithreading
    partition = get_partition(Nk.M, Threads.nthreads())

    t = 0.0
    dimension = 3
    dΩ = 4*pi

    FD = force_dissipation(Nk.M)

    interpScheeme = interp_scheeme(Nk.kk)
    time_steppingScheeme  = time_stepping_scheeme(Nk.M)

    Acoustic3D(name,Nk_arguments, ω ,Nk, Stk, F1, FSt, partition,interpScheeme,time_steppingScheeme,diags, t, c, dimension,dΩ, FD)
end

