# The structures and the constructors of the Petviashvilli runs go here

@doc raw"""
    Petviashvilli

Simulation structure for Petviashvilli wave turbulence It contains

    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``.  1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments
   
    Nk::wave_spectrum_khkz #wave action
    Sk::wave_spectrum_khkz #collisional integral
    F1::wave_spectrum_khkz #working field
    FSt::Vector{wave_spectrum} #Array of working fields for multithreading
    partition::Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}} #partition for multithreading
   
    # Type of interpolation, integration, and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    integ_scheeme::Integ_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    dimension::Int # physical dimension of the system
    dΩ::Float64 # surface of the unit sphere

    FD::force_dissipation # Contains all the terms about force and dissipation. 


"""
mutable struct Petviashvilli{Interp_Scheeme<:abstract_interpolation,Integ_Scheeme<:abstract_integration,Time_Stepping<:abstract_time_stepping}
    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``.  1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments
   
    Nk::wave_spectrum_khkz #wave action
    Sk::wave_spectrum_khkz #collisional integral
    F1::wave_spectrum_khkz #working field
    FSt::Vector{wave_spectrum} #Array of working fields for multithreading
    partition::Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}} #partition for multithreading
   
    # Type of interpolation, integration, and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    integ_scheeme::Integ_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    dimension::Int # physical dimension of the system
    dΩ::Float64 # surface of the unit sphere

    FD::force_dissipation # Contains all the terms about force and dissipation. 
end

@doc raw"""
    Petviashvilli(Nk::wave_spectrum_khkz; a=0.5, c=1.0,interp_scheeme=bilin_interp_khkz,integ_scheeme=integrate_with_log_bins_khkz,drive_scheeme=RK2_step,M=0)

Constructor of a `Petviashvilli` structure. Optionally we set interpolation and time-stepping scheemes:

    interp_scheeme : bilin_interp_khkz (default), cpow_interp_khkz
    integ_scheeme : integrate_with_log_bins_khkz (default), integrate_with_cpow_khkz
    drive_scheeme : Euler_step, RK2_step (default), ETD2_step

Note that the integration scheeme is used only for computing integrals in the `(kh, kz)` space. 
Collisional integrals are 1D and are computed using integrate_with_log_bins(). 

The parameter `a=0.5` and `c=1.0` are the dispersive length and the speed of sound

"""

function Petviashvilli(Nk::wave_spectrum_khkz;interp_scheeme=bilin_interp_khkz, integ_scheeme=integrate_with_log_bins_khkz, drive_scheeme=Euler_step,nthx=0, nthy=0)

    name = "Petviashvilli"
    Nk_arguments = 2
    ω = ω_Petviashvilli

    kkh = Nk.kkh
    khmin = kkh[1]
    khmax = kkh[end]
    Mh = Nk.Mh
    kkz = Nk.kkz
    kzmin = kkz[1]
    kzmax = kkz[end]
    Mz = Nk.Mz

    Stk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    F1 = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    FSt = Vector{wave_spectrum}(undef, Threads.nthreads())

    for ith = 1:Threads.nthreads()
        FSt[ith] = wave_spectrum(khmin,khmax,Mh)
    end

    
    diags = diagnostic_container(Mh,Mz)

    #Set partition for multithreading
    partition = get_partition(Mh,Mz,Threads.nthreads(),nthx,nthy)

    t = 0.0
    dimension = 1
    dΩ = 2.0

    FD = force_dissipation(Mh,Mz)

    interpScheeme = interp_scheeme(kkh,kkz)
    integScheeme = integ_scheeme()
    driveScheeme = drive_scheeme(Mh,Mz)

    Petviashvilli(name,Nk_arguments,ω, Nk, Stk, F1, FSt, partition, interpScheeme, 
        integScheeme, driveScheeme, diags, t, dimension, dΩ, FD)

end




@doc raw"""
    Petviashvilli_Asymp

Simulation structure for Petviashvilli for 'kx<<ky' wave turbulence. It contains

    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``.  1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments

    Nk::wave_spectrum_khkz #wave action
    Sk::wave_spectrum_khkz #collisional integral
    γk::Array{Float64,2} # gamma term of WKE
    ηk::Array{Float64,2} # gamma term of WKE

    F1::wave_spectrum_khkz #working field
    FSt::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt1::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt2::Vector{wave_spectrum} #Array of working fields for multithreading
    partition::Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}} #partition for multithreading
    
    # Type of interpolation, integration, and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    integ_scheeme::Integ_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    dimension::Int # physical dimension of the system
    dΩ::Float64 # surface of the unit sphere

    FD::force_dissipation # Contains all the terms about force and dissipation. 


"""
mutable struct Petviashvilli_Asymp{Interp_Scheeme<:abstract_interpolation,Integ_Scheeme<:abstract_integration,Time_Stepping<:abstract_time_stepping}
    name::String #name of the simulation type
    Nk_arguments::Int # Number of arguments of ``n_k``.  1: (fully symetric) , 2: (cylindrical average in 3D or mirror symmetric in 2D), 3: Only mirror symmetric in 3D
    ω # Dispersion relation. This is a function of ``k``. It takes `Nk_argument` arguments

    Nk::wave_spectrum_khkz #wave action
    Sk::wave_spectrum_khkz #collisional integral
    γk::Array{Float64,2} # gamma term of WKE
    ηk::Array{Float64,2} # gamma term of WKE

    F1::wave_spectrum_khkz #working field
    FSt::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt1::Vector{wave_spectrum} #Array of working fields for multithreading
    FSt2::Vector{wave_spectrum} #Array of working fields for multithreading
    partition::Vector{Tuple{UnitRange{Int64},UnitRange{Int64}}} #partition for multithreading
    
    # Type of interpolation, integration, and time stepping scheemes
    interp_scheeme::Interp_Scheeme
    integ_scheeme::Integ_Scheeme
    time_stepping::Time_Stepping

    # Outputs and diagnostics
    diags::diagnostic_container

    t::Float64 #current time
    dimension::Int # physical dimension of the system
    dΩ::Float64 # surface of the unit sphere

    FD::force_dissipation # Contains all the terms about force and dissipation. 
end

@doc raw"""
    Petviashvilli_Asymp(Nk::wave_spectrum_khkz; a=0.5, c=1.0,interp_scheeme=bilin_interp_khkz,integ_scheeme=integrate_with_log_bins_khkz,drive_scheeme=RK2_step,M=0)

Constructor of a `Acoustic2D` structure. Optionally we set interpolation and time-stepping scheemes:

    interp_scheeme : bilin_interp_khkz (default), cpow_interp_khkz
    integ_scheeme : integrate_with_log_bins_khkz (default), integrate_with_cpow_khkz
    drive_scheeme : Euler_step, RK2_step (default), ETD2_step

Note that the integration scheeme is used only for computing integrals in the `(kh, kz)` space. 
Collisional integrals are 1D and are computed using integrate_with_log_bins(). 

The parameter `a=0.5` and `c=1.0` are the dispersive length and the speed of sound

"""

function Petviashvilli_Asymp(Nk::wave_spectrum_khkz;interp_scheeme=bilin_interp_khkz, integ_scheeme=integrate_with_log_bins_khkz, drive_scheeme=Euler_step,nthx=-1, nthy=-1)

    name = "Petviashvilli_Asymp"
    Nk_arguments = 2
    ω = ω_Petviashvilli_Asymp

    kkh = Nk.kkh
    khmin = kkh[1]
    khmax = kkh[end]
    Mh = Nk.Mh
    kkz = Nk.kkz
    kzmin = kkz[1]
    kzmax = kkz[end]
    Mz = Nk.Mz

    Stk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    γk =  zeros(Mh,Mz)
    ηk =  zeros(Mh,Mz)

    F1 = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)
    FSt = Vector{wave_spectrum}(undef, Threads.nthreads())
    FSt1 = Vector{wave_spectrum}(undef, Threads.nthreads())
    FSt2 = Vector{wave_spectrum}(undef, Threads.nthreads())

    for ith = 1:Threads.nthreads()
        FSt[ith] = wave_spectrum(khmin,khmax,Mh)
        FSt1[ith] = wave_spectrum(khmin,khmax,Mh)
        FSt2[ith] = wave_spectrum(khmin,khmax,Mh)
    end

    diags = diagnostic_container(Mh,Mz)
    diags.glob_diag["Mx"]= WavKinS.global_ouput("Total potential enstrophy")
    diags.glob_diag["dMx"]= WavKinS.global_ouput("Total potential enstrophy dissipation")
    diags.glob_diag["Phi"]= WavKinS.global_ouput("Total zonostrophy")
    diags.glob_diag["dPhi"]= WavKinS.global_ouput("Total zonostrophy dissipation")

    #Set partition for multithreading
    partition = get_partition(Mh,Mz,Threads.nthreads(),nthx,nthy)
    
    t = 0.0
    dimension = 1
    dΩ = 2.0

    FD = force_dissipation(Mh,Mz)

    interpScheeme = interp_scheeme(kkh,kkz)
    integScheeme = integ_scheeme()
    driveScheeme = drive_scheeme(Mh,Mz)

    Petviashvilli_Asymp(name,Nk_arguments,ω,Nk, Stk,γk,ηk, F1, FSt,FSt1,FSt2, partition, interpScheeme, 
        integScheeme, driveScheeme, diags, t, dimension, dΩ, FD)

end