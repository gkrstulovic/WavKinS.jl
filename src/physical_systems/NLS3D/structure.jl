# The structures and the constructors of the NLS runs go here


@doc raw"""
    NLS3D

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
    FSt3::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt4::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt5::Vector{wave_spectrum} #Array of working fields for multithreading
    
    partition::Vector{UnitRange{Int64}} #partition for multithreading

    # Type of interpolation ,integration, and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    integ_scheeme::Integ_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    dimension::Int # physical dimension of the system
    dΩ::Float64 # surface of the unit sphere

    FD::force_dissipation # Contains all the terms about force and dissipation 

"""
mutable struct NLS3D{Interp_Scheeme<:abstract_interpolation,Integ_Scheeme<:abstract_integration,Time_Stepping<:abstract_time_stepping}
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
    FSt3::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt4::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt5::Vector{wave_spectrum} #Array of working fields for multithreading
    
    partition::Vector{UnitRange{Int64}} #partition for multithreading

    # Type of interpolation ,integration, and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    integ_scheeme::Integ_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    dimension::Int # physical dimension of the system
    dΩ::Float64 # surface of the unit sphere

    FD::force_dissipation # Contains all the terms about force and dissipation 
end

@doc raw"""
    NLS3D(Nk::wave_spectrum;interp_scheeme=lin_interp,drive_scheeme=RK2_step)

Constructor of a `NLS3D` structure. Optionally we set interpolation and time-stepping scheemes:

* `interp_scheem` : lin_interp (default), powexp_interp, powGauss_interp,BS_interp
* `time_stepping_scheeme`: [`Euler_step`](@ref), [`RK2_step`](@ref) (default), [`RK4_step`](@ref), [`ETD2_step`](@ref) and [`ETD4_step`](@ref)

The parameter `β=0.` is the exponent of non-linear term

"""
function NLS3D(Nk::wave_spectrum;interp_scheeme=lin_interp,integ_scheeme=integrate_with_log_bins_khkz,drive_scheeme=RK2_step)

    name = "NLS3D"
    Nk_arguments = 1
    ω = ω_NLS

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
    FSt3 = Vector{wave_spectrum}(undef, Threads.nthreads())
    FSt4 = Vector{wave_spectrum}(undef, Threads.nthreads())
    FSt5 = Vector{wave_spectrum}(undef, Threads.nthreads())
  
    
    for ith = 1:Threads.nthreads()
        FSt[ith] = wave_spectrum(kmin^2,kmax^2,M)
        FSt1[ith] = wave_spectrum(kmin^2,kmax^2,M)
        FSt2[ith] = wave_spectrum(kmin^2,kmax^2,M)
        FSt3[ith] = wave_spectrum(kmin^2,kmax^2,M)
        FSt4[ith] = wave_spectrum(kmin^2,kmax^2,M)
        FSt5[ith] = wave_spectrum(kmin^2,kmax^2,M)
    end

    diags = diagnostic_container(Nk.M)
    diags.sp_outs["Qk"] = spectral_output("Waveaction flux",Nk.M)
    diags.sp_store["Qk"] = Array{Float64}(undef, Nk.M, 0)


    partition = get_partition(Nk.M, Threads.nthreads())

    t = 0.0
    dimension = 3
    dΩ = 4π

    FD = force_dissipation(M)

    interpScheeme = interp_scheeme(kk)
    integScheeme = integ_scheeme()
    driveScheeme  = drive_scheeme(M)

    NLS3D(name,Nk_arguments,ω,Nk, Stk,γk,ηk, F1, FSt,FSt1,FSt2,FSt3,FSt4,FSt5, partition,interpScheeme,integScheeme,driveScheeme,
                    diags, t, dimension,dΩ, FD)
end
