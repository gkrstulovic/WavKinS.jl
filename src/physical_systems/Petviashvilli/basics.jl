# Here we define the matrices and dispersion relation and resonant manifold of Petviashvilli runs


#--- ω_Petviashvilli model
function ω_Petviashvilli(kx, ky)
    return kx * (kx^2 + ky^2)
end

function S_Petviashvilli(k1x, k2x, k3x)
    return 0.5 * abs(k1x * k2x * k3x)
end

function get_k1y_Petviashvilli(kx, ky, k1x)
    Discrim = sqrt(3 * k1x * kx^2 * (kx - k1x) + (k1x^2 - k1x * kx + kx^2) * ky^2)
    k1y_plus = -(k1x * ky - kx * ky + Discrim) / kx
    k1y_minus = -(k1x * ky - kx * ky - Discrim) / kx
    return k1y_plus, k1y_minus
end

function Det_Petviashvilli(kx, ky, k1x)
    Delta = 2 * sqrt(kx^2 * ky^2 + k1x^2 * (-3 * kx^2 + ky^2) + k1x * (3 * kx^3 - kx * ky^2))
    return Delta
end

function get_ksup_Petviashvilli(kx, ky)
    ksup = 0.5 * kx * (1.0 + sqrt(3 * (kx^2 + ky^2) / abs(ky^2 - 3 * kx^2)))
    return ksup
end

function integrable_singularity_contribution_Petviashvilli(kx, ky, kend, ksup)
    dem = 2 * (3 * kx^2 * (kx^2 + ky^2) * (3 * kx^2 - ky^2))^0.25
    return 2 * sqrt(ksup - kend) / dem
end

# ---------------------------------------------------------
#--- Petviashvilli_Asymp model
function ω_Petviashvilli_Asymp(kx, ky)
    return kx * ky^2
end

function ρ_zonostrophy_Petviashvilli_Asymp(kx, ky)
    return kx^3 / ky^2
end

function ρ_Potential_Petviashvilli_Asymp(kx, ky)
    return kx
end

function get_k1x_Petviashvilli_Asymp(kx, ky, k1y)
    return kx * k1y * (2 * ky - k1y) / (ky * (2 * k1y - ky))
end

function Det_Petviashvilli_Asymp(ky, k1y)
    Delta = abs(ky * (2 * k1y - ky) + 1.e-80)
    return Delta
end
