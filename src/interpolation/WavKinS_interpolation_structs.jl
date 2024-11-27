abstract type abstract_interpolation end

#############     Linear interpolation #############
struct lin_interp_data
    α::Vector{Float64} # Constant of linear approximation
    β::Vector{Float64} # slope of linear approximation
end

function lin_interp_data(kk::AbstractVector)
    M = length(kk)
    lin_interp_data(zeros(M), zeros(M))
end

struct lin_interp <: abstract_interpolation
    coeff::lin_interp_data
end

@doc raw"""
    lin_interp(kk::AbstractVector)

Allocate coefficients for linear interpolation and define interpolation method. See [`val_nk`](@ref).
"""
function lin_interp(kk::AbstractVector)
    lin_interp(lin_interp_data(kk))
end
########################################################

#############     Linear interpolation in log scale #############
struct linlog_interp_data
    C0::Vector{Float64}# Prefactor
    α::Vector{Float64}# Slope 
end

function linlog_interp_data(kk::AbstractVector)
    M = length(kk)
    linlog_interp_data(zeros(M), zeros(M))
end

struct linlog_interp <: abstract_interpolation
    coeff::linlog_interp_data
end

@doc raw"""
    linlog_interp(kk::AbstractVector)

Allocate coefficients for linear interpolation in log scale (power law) and define interpolation method. See [`val_nk`](@ref).
"""
function linlog_interp(kk::AbstractVector)
    linlog_interp(linlog_interp_data(kk))
end
########################################################

#############   Power-exp interpolation #############
struct powexp_interp_data # we assum n[k] = c0 k^{-α} Exp[-β k]
    c0::Vector{Float64}# Amplitud
    α::Vector{Float64}# Constant of linear approximation
    β::Vector{Float64}# slope of linear approximation
end

function powexp_interp_data(kk::AbstractVector)
    M = length(kk)
    powexp_interp_data(zeros(M), zeros(M), zeros(M))
end

struct powexp_interp <: abstract_interpolation
    coeff::powexp_interp_data
end

@doc raw"""
    powexp_interp(kk::AbstractVector)

Allocate coefficients for power-law exponential interpolation. See [`val_nk`](@ref).
"""
function powexp_interp(kk::AbstractVector)
    powexp_interp(powexp_interp_data(kk))
end
########################################################

#############   Power-exp interpolation #############
struct powGauss_interp_data # we assume n[k] = c0 k^{-α} Exp[-β k^2]
    c0::Vector{Float64}# Amplitud
    α::Vector{Float64}# Constant of linear approximation
    β::Vector{Float64}# slope of linear approximation
end

function powGauss_interp_data(kk::AbstractVector)
    M = length(kk)
    powGauss_interp_data(zeros(M), zeros(M), zeros(M))
end

struct powGauss_interp <: abstract_interpolation
    coeff::powGauss_interp_data
end

@doc raw"""
    powGauss_interp(kk::AbstractVector)

Allocate coefficients for power-law Gaussian interpolation. See [`val_nk`](@ref).
"""
function powGauss_interp(kk::AbstractVector)
    powGauss_interp(powGauss_interp_data(kk))
end
########################################################

#############  BS spline interpolation #############
using BSplineKit

mutable struct BS_interp_data{Interp}
    spline::Interp
end
function BS_interp_data(kk::AbstractVector)
    val_nk = interpolate(kk, zeros(length(kk)), BSplineOrder(4), Natural())
    BS_interp_data(val_nk)
end


struct BS_interp{Data<:BS_interp_data} <: abstract_interpolation
    coeff::Data
end




@doc raw"""
    BS_interp(kk::AbstractVector)

Allocate nodes for BSplines interpolation. See [`val_nk`](@ref).
"""
function BS_interp(kk::AbstractVector)
    BS_interp(BS_interp_data(kk))
end
#############  BS spline interpolation #############

#############  bilinear interpolation khkz #############
struct bilin_interp_data_khkz # we assume n(kh, kz) = c0 - αh kh - αz kz - β kh 
    c0::Array{Float64,2} # constant
    αh::Array{Float64,2} # -slope along horizontal 
    αz::Array{Float64,2} # -slope along vertical 
    β::Array{Float64,2} # nonlinear term coefficient
end

function bilin_interp_data_khkz(kkh::AbstractVector, kkz::AbstractVector)
    Mh = length(kkh)
    Mz = length(kkz)
    bilin_interp_data_khkz(zeros(Mh, Mz), zeros(Mh, Mz), zeros(Mh, Mz), zeros(Mh, Mz))
end

struct bilin_interp_khkz <: abstract_interpolation
    coeff::bilin_interp_data_khkz
end

@doc raw"""
    bilin_interp_khkz(kkh::AbstractVector, kkz::AbstractVector)

Allocate coefficients for bilinear interpolation and define interpolation method. See [`val_nk`](@ref).
"""
function bilin_interp_khkz(kkh::AbstractVector, kkz::AbstractVector)
    bilin_interp_khkz(bilin_interp_data_khkz(kkh, kkz))
end
########################################################

#############  cpow interpolation khkz #############
struct cpow_interp_data_khkz # we assume n(kh, kz) = β + c0 kh^(-αh) kz^(-αz) 
    c0::Array{Float64,2}
    αh::Array{Float64,2}
    αz::Array{Float64,2}
    β::Array{Float64,2}
end

function cpow_interp_data_khkz(kkh::AbstractVector, kkz::AbstractVector)
    Mh = length(kkh)
    Mz = length(kkz)
    cpow_interp_data_khkz(zeros(Mh, Mz), zeros(Mh, Mz), zeros(Mh, Mz), zeros(Mh, Mz))
end

struct cpow_interp_khkz <: abstract_interpolation
    coeff::cpow_interp_data_khkz
end

@doc raw"""
    cpow_interp_khkz(kkh::AbstractVector, kkz::AbstractVector)

Allocate coefficients for cpow interpolation and define interpolation method. See [`val_nk`](@ref).
"""
function cpow_interp_khkz(kkh::AbstractVector, kkz::AbstractVector)
    cpow_interp_khkz(cpow_interp_data_khkz(kkh, kkz))
end
########################################################

#############  bilinear interpolation in log khkz #############
struct bilinlog_interp_data_khkz # we assume log (n(kh, kz)) = c0 - αh log(kh) - αz log(kz) - β log(kh) log(kz) 
    c0::Array{Float64,2} # constant in log scale
    αh::Array{Float64,2} # -slope along horizontal in log scale
    αz::Array{Float64,2} # -slope along vertical in log scale
    β::Array{Float64,2} # nonlinear term coefficient in log scale
end

function bilinlog_interp_data_khkz(kkh::AbstractVector, kkz::AbstractVector)
    Mh = length(kkh)
    Mz = length(kkz)
    bilinlog_interp_data_khkz(zeros(Mh, Mz), zeros(Mh, Mz), zeros(Mh, Mz), zeros(Mh, Mz))
end

struct bilinlog_interp_khkz <: abstract_interpolation
    coeff::bilinlog_interp_data_khkz
end

@doc raw"""
    bilinlog_interp_khkz(kkh::AbstractVector, kkz::AbstractVector)

Allocate coefficients for bilinear interpolation in log scale and define interpolation method. See [`val_nk`](@ref).
"""
function bilinlog_interp_khkz(kkh::AbstractVector, kkz::AbstractVector)
    bilinlog_interp_khkz(bilinlog_interp_data_khkz(kkh, kkz))
end
########################################################

#############  exponential interpolation khkz #############
struct exp_interp_data_khkz # we assume log (n(kh, kz)) = c0 - αh kh - αz kz - β kh kz
    c0::Array{Float64,2} # constant in log scale
    αh::Array{Float64,2} # -slope along horizontal in log scale
    αz::Array{Float64,2} # -slope along vertical in log scale
    β::Array{Float64,2} # nonlinear term coefficient in log scale
end

function exp_interp_data_khkz(kkh::AbstractVector, kkz::AbstractVector)
    Mh = length(kkh)
    Mz = length(kkz)
    exp_interp_data_khkz(zeros(Mh, Mz), zeros(Mh, Mz), zeros(Mh, Mz), zeros(Mh, Mz))
end

struct exp_interp_khkz <: abstract_interpolation
    coeff::exp_interp_data_khkz
end

@doc raw"""
    exp_interp_khkz(kkh::AbstractVector, kkz::AbstractVector)

Allocate coefficients for exponential interpolation and define interpolation method. See [`val_nk`](@ref).
"""
function exp_interp_khkz(kkh::AbstractVector, kkz::AbstractVector)
    exp_interp_khkz(exp_interp_data_khkz(kkh, kkz))
end
########################################################