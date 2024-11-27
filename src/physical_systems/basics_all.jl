# Here we define quantities that can be useful for all runs

@doc raw"""
    Delta(k::Float64, k1::Float64, k2::Float64)

Assuming that ``k_1 < k + k_2``, ``k_2 < k + k_1``, and ``k < k_1 + k_2``, it returns twice the area of the triangle of sides ``k``, ``k_1``, and ``k_2``, i.e. ``\Delta = \frac{1}{2} \sqrt{(k + k_1 + k_2) (k + k_1 - k_2) (k - k_1 + k_2) (-k + k_1 + k_2)} = 2 \left[ (k k_1)^2 + (k k_2)^2 + (k_1 k_2)^2 \right] - k^4 - k_1^4 - k_2^4``. 
"""
function Delta(k::Float64, k1::Float64, k2::Float64)
    x = 2 * ((k * k1)^2 + (k * k2)^2 + (k1 * k2)^2) - k^4 - k1^4 - k2^4
    return 0.5 * sqrt(x)
end

@doc raw"""
    Δ_vs(k::Float64, p::Float64, q::Float64)

Return ``\Delta = \frac{1}{2} \sqrt{(k^2-p^2) q (2 k+q)}``. It is [`Delta`](@ref) with ``k_1 = (k + p + q)/2`` and ``k_2 = (k - p + q)/2``.
"""
function Δ_vs(k::Float64, p::Float64, q::Float64)
    return sqrt((k^2-p^2)*q*(2*k+q))/2
end


@doc raw"""
    θkmax(k,kmax)

Gives ones if ``k < k_{max}``, ``0``, otherwise. 
"""
function θkmax(k,kmax)
    if  abs(k)<=kmax
        return 1.
    else
        return 0.
    end
end