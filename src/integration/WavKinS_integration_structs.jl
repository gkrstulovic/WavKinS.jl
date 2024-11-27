abstract type abstract_integration end


@doc raw"""
    integrate_with_log_bins_khkz()

Select the use of the trapezoidal integration with log bins khkz method. See [`integrate`](@ref).

"""
struct integrate_with_log_bins_khkz<:abstract_integration end


@doc raw"""
    integrate_with_cpow_khkz()

Select the use of the cpow integration with log bins khkz method. See [`integrate`](@ref).

"""
struct integrate_with_cpow_khkz<:abstract_integration end



