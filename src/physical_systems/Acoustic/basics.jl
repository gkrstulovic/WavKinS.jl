# Here we define the matrices and dispersion relation and resonant manifold of Acoustic runs


# -------- 2D isotropic acoustic turbulence: --------------
function S_2DAcoustic(k1, k2, k3; V0=1.0)
    return V0^2 * k1 * k2 * k3
end

function S_2DAcoustic_over_k(k2, k3; V0=1.0)
    return V0^2 * k2 * k3
end

function ω_Acoustic(k; c=1.0)
    return c * k
end

function ω_Acoustic_2args(kx,ky; c=1.0)
    return ω_Acoustic(sqrt(kx^2+ky^2); c)
end

function ω_Acoustic_disp(k; c=1.0, a=0.5)
    return c * k * (1 + (a * k)^2)
end

# ---------------------------------------------------------

# -------- 3D isotropic acoustic turbulence: --------------
function S_3DAcoustic(k1, k2, k3; V0=1.0)
    return V0^2 * k1^2 * k2^2 * k3^2
end

function S_3DAcoustic_over_k(k2, k3; V0=1.0)
    return V0^2 * k2^2 * k3^2
end