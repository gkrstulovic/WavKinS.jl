@doc raw"""
    cardano(p::Float64, q::Float64)

Assuming that `p` > 0 and `q` < 0, it returns the positive real root of ``x^3 + `` `p` `` x + `` `q`. 

"""
function cardano(p::Float64, q::Float64)
    d = (q^2/4.0 + p^3/27.0)^(1/2)
    x = (-0.5*q + d)^(1/3) - (0.5*q + d)^(1/3)
    return max(x,0.0)
end

@doc raw"""
    NewtonRaphson(f, x0; ϵ=1e-9, tol=1e-12, maxIter = 1000)

Find the zero of the function `f` using the Newton-Raphson algorithm. 

* `x0`: starting point
* `ϵ`: spacing used to estimate derivative
* `tol`: tolerence
* `maxIter`: maximal number of iterations

Adapted from [here](https://www.matecdev.com/posts/julia-newton-raphson.html#newton-raphson-method).

"""
function NewtonRaphson(f, x0; ϵ=1e-9, tol=1e-12, maxIter = 1000)
    x = x0
    fx = f(x0)
    iter = 0
    while abs(fx) > tol && iter < maxIter
        fpx = (f(x+ϵ/2) - f(x-ϵ/2)) / ϵ # numerical estimate of the derivative 
        x = x  - fx/fpx     
        fx = f(x)          
        iter += 1
    end
    return x
end

@doc raw"""
    NewtonRaphson_with_derivative(f, fp, x0, tol=1e-12, maxIter = 1000)

Same as [`NewtonRaphson`](@ref "NewtonRaphson"), but using the explicit derivative of `f`, which is the argument `fp`, instead of a numerical estimate.

"""
function NewtonRaphson_with_derivative(f, fp, x0; tol=1e-12, maxIter = 1000)
    x = x0
    fx = f(x0)
    iter = 0
    while abs(fx) > tol && iter < maxIter
        fpx = fp(x)         # explicit evaluation of the derivative 
        x = x  - fx/fpx     
        fx = f(x)           
        iter += 1
    end
    return x
end

