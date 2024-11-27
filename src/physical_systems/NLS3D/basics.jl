# Here we define the matrices and dispersion relation and resonant manifold of NLS 3D runs

function ω_NLS(k)
    return abs(k)^2
end

function T1234squared_NLS3D(k1, k2, k3, k4)
    return min(k1, k2, k3, k4) * k1 * k2 * k3 * k4

end

function T1234squared_NLS3D_in_ω(ω1, ω2, ω3, ω4)
    return sqrt(min(ω1, ω2, ω3, ω4))

end


function get_k1_NLS(k, k2, k3)
    return sqrt(k2^2 + k3^2 - k^2)
end

function Det_NLS3D(k, k1, k2, k3)
    return 2 * k1
end
