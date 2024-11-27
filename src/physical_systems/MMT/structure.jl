# The structures and the constructors of the MMT runs go here


@doc raw"""
    MMT

Simulation structure for MMT wave turbulence. It contains

    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``.  1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments
    
    Nk::wave_spectrum #wave action
    Sk::wave_spectrum #collisional integral
    γk::Vector{Float64} # gamma term of WKE
    ηk::Vector{Float64} # eta term of WKE

    F1::wave_spectrum #working field
    FSt::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt1::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt2::Vector{wave_spectrum} #Array of working fields for multithreading
    partition::Vector{UnitRange{Int64}} #partition for multithreading

    # Type of interpolation and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    β::Float64 # dispersive length
    dimension::Int # physical dimension of the system
    dΩ::Float64 # surface of the unit sphere

    FD::force_dissipation # Contains all the terms about force and dissipation 


"""
mutable struct MMT{Interp_Scheeme<:abstract_interpolation,Time_Stepping<:abstract_time_stepping}
    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``.  1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments
    
    Nk::wave_spectrum #wave action
    Sk::wave_spectrum #collisional integral
    γk::Vector{Float64} # gamma term of WKE
    ηk::Vector{Float64} # eta term of WKE

    F1::wave_spectrum #working field
    FSt::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt1::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt2::Vector{wave_spectrum} #Array of working fields for multithreading
    partition::Vector{UnitRange{Int64}} #partition for multithreading

    # Type of interpolation and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    β::Float64 # dispersive length
    dimension::Int # physical dimension of the system
    dΩ::Float64 # surface of the unit sphere

    FD::force_dissipation # Contains all the terms about force and dissipation 
end

@doc raw"""
    MMT(Nk::wave_spectrum; β=0.,interp_scheeme=lin_interp,drive_scheeme=RK2_step

Constructor of a `MMT` structure. Optionally we set interpolation and time-stepping scheemes:

    interp_scheeme : lin_interp (default), powexp_interp, powGauss_interp,BS_interp
    drive_scheeme : Euler_step, RK2_step (default), ETD2_step, AB_Euler_step, and AB2_RK2_step.

The parameter `β=0.` is the exponent of non-linear term

"""
function MMT(Nk::wave_spectrum; β=0.,interp_scheeme=lin_interp,drive_scheeme=RK2_step)

    name = "MMT"
    Nk_arguments = 1
    ω = ω_MMT

    kk = Nk.kk
    kmin = kk[1]
    kmax = kk[end]
    M = Nk.M

    Stk = wave_spectrum(kmin,kmax,M)
    γk =  zeros(M)
    ηk =  zeros(M)

    F1 = wave_spectrum(kmin,kmax,M)
    FSt = Vector{wave_spectrum}(undef, Threads.nthreads())
    FSt1 = Vector{wave_spectrum}(undef, Threads.nthreads())
    FSt2 = Vector{wave_spectrum}(undef, Threads.nthreads())
    for ith = 1:Threads.nthreads()
        FSt[ith] = wave_spectrum(kmin,kmax,M)
        FSt1[ith] = wave_spectrum(kmin,kmax,M)
        FSt2[ith] = wave_spectrum(kmin,kmax,M)
    end

    diags = diagnostic_container(Nk.M)
    diags.sp_outs["Qk"] = spectral_output("Waveaction flux",Nk.M)
    diags.sp_store["Qk"] = Array{Float64}(undef, Nk.M, 0)


    partition = get_partition(Nk.M, Threads.nthreads())

    t = 0.0
    dimension = 1
    dΩ = 2.

    FD = force_dissipation(M)

    interpScheeme = interp_scheeme(kk)
    driveScheeme  = drive_scheeme(M)

    MMT(name,Nk_arguments,ω,Nk, Stk,γk,ηk, F1, FSt,FSt1,FSt2, partition,interpScheeme,driveScheeme,
                    diags, t, β, dimension,dΩ, FD)
end