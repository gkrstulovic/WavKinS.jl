# Here we define the matrices and dispersion relation and resonant manifold of MMT runs
using Trapz

function ω_MMT(k)
    return sqrt(abs(k))
end

function T1234squared_MMT(k1, k2, k3, k4; β=0.0)
    return abs(k1 * k2 * k3 * k4)^(β / 2)

end

function get_k1_MMT(k, k3)
    q = k3 / k
    if q < -1
        return -k * (-1 + sqrt(-q) + q)^2 / (-1 + sqrt(-q))^2
    elseif -1 <= q <= 0
        return k * q * (-1 - 2 * sqrt(-q) + q) / (1 + q)^2
    elseif 0 <= q <= 1
        return 0.5 * k * (-1 + q + sqrt(1 - 6 * q + 8 * q^(3 / 2) - 3 * q^2))
    elseif 1 <= q
        return 0.5 * k * (-1 + q + sqrt(-3 + 8 * sqrt(q) - 6 * q + q^2))
    else
        return NaN
    end
end

function get_k3_MMT(k, k1)
    q = k1 / k
    if q < -1
        return -k * (-1 + sqrt(-q) + q)^2 / (-1 + sqrt(-q))^2
    elseif -1 <= q <= 0
        return k * q * (-1 - 2 * sqrt(-q) + q) / (1 + q)^2
    elseif 0 <= q <= 1
        return 0.5 * k * (-1 + q + sqrt(1 - 6 * q + 8 * q^(3 / 2) - 3 * q^2))
    elseif 1 <= q
        return 0.5 * k * (-1 + q + sqrt(-3 + 8 * sqrt(q) - 6 * q + q^2))
    else
        return NaN
    end
end

function Det_MMT(k,k1, k2, k3)
    eps = 1e-80
    Delta = -sign(k2) / sqrt(abs(k2) + eps) + sign(k1) / sqrt(abs(k1) + eps)
    return abs(Delta) / 2.0
end



## Computation of the derivative of the dimensionless collisional integral. 

function dxdIgen(x, q1, q2, q3; β=0.0)
    dxI = log(abs(q1 * q2 * q3)) * abs(q1 * q2 * q3)^(0.5 * β - x) * (abs(q1)^x + abs(q2)^x - abs(q3)^x - 1.0) / WavKinS.Det_MMT(q1, q2, q3, 1.0)
    return dxI - abs(q1 * q2 * q3)^(0.5 * β - x) * (log(abs(q1)) * abs(q1)^x + log(abs(q2)) * abs(q2)^x - log(abs(q3)) * abs(q3)^x) / WavKinS.Det_MMT(q1, q2, q3, 1.0)
end

function q3(q1)
    return get_k3_MMT(1.0, q1)
end

function dxdI(x, q1; β=0.0)
    return dxdIgen(x, q1, q3(q1) + 1.0 - q1, q3(q1); β=β)
end


@doc raw"""
    dxI_MMT(x; β=0.0, kmax=10000, npoints=100000)

Compute the derivative of the dimensionless colissional integral for the MMT model
*  `x`: value of the exponent (argument of the funcion)
* `β` : homogeneity degree of the interaction matrix. (default `β=0`) 
* `kmax` : Ultraviolet cut off. (default `kmax=10000`) 
* `npoints` : number of discretisation points. (default `npoints=100000`) 
"""
function dxI_MMT(x; β=0.0, kmax=10000, npoints=100000)

    Q1 = range(-kmax, kmax, npoints)
    dIvec = similar(Q1)
    for (ind, q1) in enumerate(Q1)
        dIvec[ind] = dxdI(x, q1; β=β)
    end
    return trapz(Q1, dIvec)
end