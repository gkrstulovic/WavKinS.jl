# Here we define the matrices and dispersion relation and resonant manifold for Bogoliubov wave turbulence


function ω_Bogoliubov(k; c=1.0, ξ=1.0)
    return c * k * sqrt(1 + 0.5 * (ξ * k)^2)
end

function disp_corr_bogo(k; ξ=1.0)
    return sqrt(1 + 0.5 * (ξ * k)^2)
end


# ---------------------------------------------------------

# -------- 3D isotropic Bogoliubov turbulence: --------------
function V123squared_Bogo3D(k1, k2, k3; V0=1.0, ξ=1.0)
    γ1 = disp_corr_bogo(k1; ξ=ξ)
    γ2 = disp_corr_bogo(k2; ξ=ξ)
    γ3 = disp_corr_bogo(k3; ξ=ξ)
    V123squared = V0^2 * k1^2 * k2^2 * k3^2
    V123squared *= (0.5 + ((γ1 * γ2 * γ3) / (6 * k1 * k2 * k3 + 1e-90)) * (k1^3 / γ1 - k2^3 / γ2 - k3^3 / γ3))^2
    return V123squared / (γ1 * γ2 * γ3)  
end


# -- Resonant manifold ---
function k2_Bogoliubov(k, k1; ξ=1.0)
    ωk = ω_Bogoliubov(k; c=1.0, ξ=ξ)
    ωk1 = ω_Bogoliubov(k1; c=1.0, ξ=ξ)
    return sqrt(sqrt(1.0 + 2 * ξ^2 * (ωk - ωk1)^2) - 1.0) / ξ
end


function Det_Bogo(k2; c=1, ξ=1)
    return c * (1 + (k2 * ξ)^2) / sqrt(1 + 0.5 * (k2 * ξ)^2)
end
