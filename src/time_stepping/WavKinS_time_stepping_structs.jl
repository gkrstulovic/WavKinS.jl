abstract type abstract_time_stepping end


@doc raw"""
    Euler_step()

Set the use of a simple Euler. See [`advance!`](@ref).
"""
struct Euler_step<:abstract_time_stepping end 

function Euler_step(M)
    Euler_step()
end
function Euler_step(Mh,Mz)
    Euler_step()
end

@doc raw"""
    AB_Euler_step()

Set the use of split Exponential-Adams-Bashforth and Euler. See [`advance!`](@ref).
"""
struct AB_Euler_step<:abstract_time_stepping end 

function AB_Euler_step(M)
    AB_Euler_step()
end
function AB_Euler_step(Mh,Mz)
    AB_Euler_step()
end
####################################################################################################################################

#########################################################################################################################################
#####            split Second-order-Exponential-Adams-Bashforth and RK2  #############
@doc raw"""
    AB2_RK2arrays

Structure for split Second-order-Exponential-Adams-Bashforth and Runge-Kutta 2
"""
struct AB2_RK2arrays
    γk::VecOrMat{Float64}  #coefficients for linear term
    γkold::VecOrMat{Float64}  #coefficients for linear term
    nAB2half::VecOrMat{Float64} #temporal array for half of the evolution
end
function AB2_RK2arrays(M::Int)
    AB2_RK2arrays(zeros(M),zeros(M),zeros(M))
end
function AB2_RK2arrays(Mh::Int,Mz::Int)
    AB2_RK2arrays(zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz))
end


@doc raw"""
    AB2_RK2_step()

Set the use of split Second-order-Exponential-Adams-Bashforth and Runge-Kutta 2. See [`advance!`](@ref).
"""
struct AB2_RK2_step<:abstract_time_stepping 
    arrays::AB2_RK2arrays
end 

function AB2_RK2_step(M::Int)
    arrays = AB2_RK2arrays(M)
    AB2_RK2_step(arrays)
end
function AB2_RK2_step(Mh::Int,Mz::Int)
    arrays = AB2_RK2arrays(Mh,Mz)
    AB2_RK2_step(arrays)
end
####################################################################################################################################

#############     Runge-Kutta 2 #############
@doc raw"""
    RK2arrays

Structure for standard Runge-Kutta 2
"""
struct RK2arrays
    Nkrkh::VecOrMat{Float64}  #coefficients for linear term
end
function RK2arrays(M::Int)
    RK2arrays(zeros(M))
end
function RK2arrays(Mh::Int,Mz::Int)
    RK2arrays(zeros(Mh,Mz))
end


@doc raw"""
    RK2_step()

Set the use of a Runge-Kutta 2. See [`advance!`](@ref).
"""
struct RK2_step<:abstract_time_stepping 
    arrays::RK2arrays
end 

function RK2_step(M::Int)
    arrays = RK2arrays(M)
    RK2_step(arrays)
end
function RK2_step(Mh::Int,Mz::Int)
    arrays = RK2arrays(Mh,Mz)
    RK2_step(arrays)
end
########################################################
#
#############    Runge-Kutta 4 #############
@doc raw"""
    RK4arrays

Structure for standard Runge-Kutta 4, no forcing nor dissipation
"""
struct RK4arrays
    nk0::VecOrMat{Float64} #temporal arrays for RK scheme
    Sk1::VecOrMat{Float64} #temporal arrays for RK scheme
    Sk2::VecOrMat{Float64} #temporal arrays for RK scheme
    Sk3::VecOrMat{Float64} #temporal arrays for RK scheme
    Sk4::VecOrMat{Float64} #temporal arrays for RK scheme
end
function RK4arrays(M::Int)
    RK4arrays(zeros(M),zeros(M),zeros(M),zeros(M), zeros(M))
end
function RK4arrays(Mh::Int,Mz::Int)
    RK4arrays(zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz))
end

@doc raw"""
    RK4_step()

Set the use of 4th-order Runge-Kutta 4, no forcing nor dissipation. See [`advance!`](@ref).
"""
struct RK4_step<:abstract_time_stepping 
    arrays::RK4arrays
end

function RK4_step(M::Int)
    arrays = RK4arrays(M)
    RK4_step(arrays)
end
function RK4_step(Mh::Int,Mz::Int)
    arrays = RK4arrays(Mh,Mz)
    RK4_step(arrays)
end

#############     ETD2 #############
@doc raw"""
    ETD2arrays

Structure for ETD2 
"""
struct ETD2arrays
    Sktmp::VecOrMat{Float64}  #temporary array
    lk::VecOrMat{Float64}  #coefficients for linear term
    cnl1::VecOrMat{Float64} #coefficients for nonlinear term
    cnl2::VecOrMat{Float64} #coefficients for nonlinear term
end
function ETD2arrays(M::Int)
    ETD2arrays(zeros(M),zeros(M),zeros(M),zeros(M))
end
function ETD2arrays(Mh::Int,Mz::Int)
    ETD2arrays(zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz))
end

@doc raw"""
    ETD2_step()

Set the use of ETD2. See [`advance!`](@ref).
"""
struct ETD2_step<:abstract_time_stepping 
    arrays::ETD2arrays
end
function ETD2_step(M::Int)
    arrays = ETD2arrays(M)
    ETD2_step(arrays)
end
function ETD2_step(Mh::Int,Mz::Int)
    arrays = ETD2arrays(Mh,Mz)
    ETD2_step(arrays)
end

#############     ETD4 #############
@doc raw"""
    ETD4arrays

Structure for ETD4 
"""
struct ETD4arrays
    nk0::VecOrMat{Float64}  #coefficients for linear term
    lk1::VecOrMat{Float64}  #coefficients for linear term
    lk2::VecOrMat{Float64}  #coefficients for linear term
    cnl1::VecOrMat{Float64} #coefficients for nonlinear term
    cnl2::VecOrMat{Float64} #coefficients for nonlinear term
    cnl3::VecOrMat{Float64} #coefficients for nonlinear term
    cnl4::VecOrMat{Float64} #coefficients for nonlinear term
    Sk1::VecOrMat{Float64} #temporal arrays for RK scheme
    Sk2::VecOrMat{Float64} #temporal arrays for RK scheme
    Sk3::VecOrMat{Float64} #temporal arrays for RK scheme
    Sk4::VecOrMat{Float64} #temporal arrays for RK scheme
    fs1::VecOrMat{Float64} #temporal arrays for RK scheme
end
function ETD4arrays(M::Int)
    ETD4arrays(zeros(M),zeros(M),zeros(M),zeros(M), zeros(M),zeros(M),zeros(M),zeros(M),zeros(M),zeros(M),zeros(M),zeros(M))
end
function ETD4arrays(Mh::Int,Mz::Int)
    ETD4arrays(zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz), zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz),zeros(Mh,Mz))
end

@doc raw"""
    ETD4_step()

Set the use of ETD4. See [`advance!`](@ref).
"""
struct ETD4_step<:abstract_time_stepping 
    arrays::ETD4arrays
end

function ETD4_step(M::Int)
    arrays = ETD4arrays(M)
    ETD4_step(arrays)
end
function ETD4_step(Mh::Int,Mz::Int)
    arrays = ETD4arrays(Mh,Mz)
    ETD4_step(arrays)
end
