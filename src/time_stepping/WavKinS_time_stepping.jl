@doc raw"""
    init_temporal_scheeme!(::Euler_step, Run, dt)

Initialization of the temporal scheeme.
* First argument is the time-stepping method: `Euler_step`, `AB_Euler_step`, `AB2_RK2_step`, `RK2_step`, `ETD2_step`, `RK4_step`, `ETD4_step`.
* `Run`: run structure
* `dt`: time step
"""
function init_temporal_scheeme!(::Euler_step, Run, dt)
    println("Using Euler temporal scheeme")
    println("")
end

@doc raw"""
    advance!(::Euler_step, Run, dt)

Make one time step with `Euler_step`.
* `Run`: run structure
* `dt`: time step

The time advancement is simply
``n_{\bf k}(t+\mathrm{d}t) = \left[ St_{\rm k} + f_{\bf k} - d_{\bf k} n_{\bf k} \right] ~ \mathrm{d}t``.
"""
function advance!(::Euler_step, Run, dt)
    nk = Run.Nk.nk
    Sk = Run.Sk.nk
    f = Run.FD.f
    D = Run.FD.D

    St_k!(Run)
    @. nk = nk + dt * Sk
    @. nk = nk + dt * f - dt * D * nk
    clean_waveaction!(Run.Nk)
    Run.t += dt
    return nothing
end

function init_temporal_scheeme!(::AB_Euler_step, Run, dt)
    println("Using Split Exponential-Adams-Bashforth and Euler temporal scheeme")
    println("")
end

@doc raw"""
    advance!(::AB_Euler_step, Run, dt)

Make one time step with `AB_Euler_step`.
* `Run`: run structure
* `dt`: time step

The collision integral is decomposed as ``St_{\bf k} = - \gamma_{\bf k} n_{\bf k} + \eta_{\bf k}``. Time advancement is done with a splitting method, treating the term ``-(\gamma_{\bf k} + d_{\bf k}) n_{\bf k}`` implicitly and the term ``\eta_{\bf k} + f_{\bf k}`` explicitly with Euler method.

Note: The `Run` structure must have the coefficients ``\gamma_{\bf k}`` and ``\eta_{\bf k}`` that must be computed in `St_k!`. This scheeme is not implemented for all physical systems.
"""
function advance!(::AB_Euler_step, Run, dt)
    nk = Run.Nk.nk
    γk = Run.γk
    ηk = Run.ηk
    f = Run.FD.f
    D = Run.FD.D

    St_k!(Run; compute_Sk=false, compute_γk=true, compute_ηk=false)
    @. nk = exp(-dt * (γk + D)) * nk
    St_k!(Run; compute_Sk=true, compute_γk=false, compute_ηk=true)
    @. nk = nk + dt * (ηk + f)
    clean_waveaction!(Run.Nk)
    Run.t += dt
    return nothing
end

function init_temporal_scheeme!(temp_scheeme::AB2_RK2_step, Run, dt)
    println("Using Split Exponential-Adams-Bashforth-2 and RK2 temporal scheeme")
    println("WARNING: maybe still bugged")
    println("")
    nk = Run.Nk.nk
    γk = Run.γk
    ηk = Run.ηk
    f = Run.FD.f
    D = Run.FD.D

    St_k!(Run; compute_Sk=false, compute_γk=true, compute_ηk=false)
    temp_scheeme.arrays.γkold .= Run.γk
    @. nk = exp(-dt * (γk + D)) * nk
    St_k!(Run; compute_Sk=true, compute_γk=false, compute_ηk=true)
    @. nk = nk + dt * (ηk + f)
    Run.t += dt
end

@doc raw"""
    advance!(::AB2_RK2_step, Run, dt)

Make one time step with `AB2_RK2_step`.
* `Run`: run structure
* `dt`: time step

The collision integral is decomposed as ``St_{\bf k} = - \gamma_{\bf k} n_{\bf k} + \eta_{\bf k}``. Time advancement is done with a splitting method, treating the term ``-(\gamma_{\bf k} + d_{\bf k}) n_{\bf k}`` implicitly and the term ``\eta_{\bf k} + f_{\bf k}`` explicitly with Runge-Kutta 2 method.

Note: The `Run` structure must have the coefficients ``\gamma_{\bf k}`` and ``\eta_{\bf k}`` that must be computed in `St_k!`. This scheeme is not implemented for all physical systems.
"""
function advance!(temp_scheeme::AB2_RK2_step, Run, dt)
    nk = Run.Nk.nk
    ηk = Run.ηk
    f = Run.FD.f
    D = Run.FD.D
    γk = temp_scheeme.arrays.γk
    γkold = temp_scheeme.arrays.γkold
    nAB2half = temp_scheeme.arrays.nAB2half

    St_k!(Run; compute_Sk=false, compute_γk=true, compute_ηk=true)
    γk .= Run.γk
    @. nAB2half = exp(-dt * ((5.0 / 8.0) * γk - (1.0 / 8.0) * γkold + 0.5 * D)) * nk # first half step linear part with AB
    # @. nAB2half = nAB2half + dt * ( ηk + f) # half step of RK2 of non-linear part

    @. nk = nAB2half + 0.5 * dt * (ηk + f) # half step of RK2 of non-linear part
    St_k!(Run; compute_Sk=true, compute_γk=true, compute_ηk=true)

    @. nAB2half = nAB2half + dt * (ηk + f)  #half step of RK2 of non-linear part
    @. nk = exp(-dt * ((7.0 / 8.0) * γk - (3.0 / 8.0) * γkold + 0.5 * D)) * nAB2half # final half step linear part with AB
    #@. nk = exp( -dt * (-(11. / 24.) * γk + (5. / 72.) * γkold  + (8. / 9.) * Run.γk + .5 * D)) * nAB2half # final half step linear part with AB

    @. γkold = γk

    Run.t += dt
    return nothing
end


function init_temporal_scheeme!(::RK2_step, Run, dt)
    println("Using implicit-explicit RK2 temporal scheeme")
    println("")
end

@doc raw"""
    advance!(::RK2_step, Run, dt)

Make one time step with `RK2_step`.
* `Run`: run structure
* `dt`: time step

The time advancement is done with Runge-Kutta 2 method.
"""
function advance!(temp_scheeme::RK2_step, Run, dt)
    nk = Run.Nk.nk
    nktmp = temp_scheeme.arrays.Nkrkh
    Sk = Run.Sk.nk
    f = Run.FD.f
    D = Run.FD.D


    @. nk *= (1.0 - 0.25 * dt * D) / (1.0 + 0.25 * dt * D) # first half implicit time step (Trapezoidal rule )
    #@. nk *= exp(- 0.5 * dt * D)  # first half implicit time step (Trapezoidal rule )


    St_k!(Run)  # explicit RK2 
    @. nktmp = nk
    @. nk = nk + 0.5 * dt * (Sk + f)
    St_k!(Run)
    @. nk = nktmp + dt * (Sk + f)

    @. nk *= (1.0 - 0.25 * dt * D) / (1.0 + 0.25 * dt * D) # second half implicit time step (Trapezoidal rule )
    #@. nk *= exp(- 0.5 * dt * D) # second half implicit time step (Trapezoidal rule )

    #clean_waveaction!(Run.Nk)

    Run.t += dt
end



function init_temporal_scheeme!(temp_scheeme::ETD2_step, Run, dt)
    println("Using ETD2 temporal scheeme")
    println("")
    lk = temp_scheeme.arrays.lk
    cnl1 = temp_scheeme.arrays.cnl1
    cnl2 = temp_scheeme.arrays.cnl2

    kk = Run.Nk.kk
    D = Run.FD.D

    cc = -D
    ch = dt * cc
    c2h = dt * cc .^ 2
    @. lk = exp(ch)
    #    @. cnl1 = expm1(ch) / cc
    nl1 = dt * map(x -> abs(x) > 5.e-2 ? expm1(x) / x : 1 + x / 2 + x^2 / 6 + x^3 / 24 + x^4 / 120 + x^5 / 720 + x^6 / 5040 + x^7 / 40320 + x^8 / 362880, ch)
    nl2 = dt * map(x -> abs(x) > 5.e-2 ? (expm1(x) - x) / x^2 : 1 / 2 + x / 6 + x^2 / 24 + x^3 / 120 + x^4 / 720 + x^5 / 5040 + x^6 / 40320 + x^7 / 362880, ch)
    @. cnl1 = nl1
    @. cnl2 = nl2
    return nothing
end

@doc raw"""
    advance!(::ETD2_step, Run, dt)

Make one time step with `ETD2_step`.
* `Run`: run structure
* `dt`: time step

"""
function advance!(temp_scheeme::ETD2_step, Run, dt)
    nk = Run.Nk.nk
    Sktmp = temp_scheeme.arrays.Sktmp
    Sk = Run.Sk.nk
    f = Run.FD.f
    lk = temp_scheeme.arrays.lk
    cnl1 = temp_scheeme.arrays.cnl1
    cnl2 = temp_scheeme.arrays.cnl2


    St_k!(Run)
    @. Sktmp = Run.Sk.nk
    @. nk = lk * nk + cnl1 * (Sk + f)

    St_k!(Run)
    @. nk = nk + cnl2 * (Sk - Sktmp)

    # clean_waveaction!(Run.Nk)
    Run.t += dt
    return nothing
end

function init_temporal_scheeme!(temp_scheeme::RK4_step, Run, dt)
    println("Using RK4  temporal scheeme")
    println("")

    return nothing
end

@doc raw"""
    advance!(::RK4_step, Run, dt)

Make one time step with `RK4_step`.
* `Run`: run structure
* `dt`: time step

The time advancement is done with Runge-Kutta 4 method.
"""
function advance!(temp_scheeme::RK4_step, Run, dt)
    nk = Run.Nk.nk
    f = Run.FD.f
    D = Run.FD.D

    nk0 = temp_scheeme.arrays.nk0
    Sk1 = temp_scheeme.arrays.Sk1
    Sk2 = temp_scheeme.arrays.Sk2
    Sk3 = temp_scheeme.arrays.Sk3
    Sk4 = temp_scheeme.arrays.Sk4

    @. nk0 = nk

    St_k!(Run)
    @. Sk1 = Run.Sk.nk + f - D * nk
    @. nk = nk0 + 0.5 * dt * Sk1


    St_k!(Run)
    @. Sk2 = Run.Sk.nk + f - D * nk
    @. nk = nk0 + 0.5 * dt * Sk2


    St_k!(Run)
    @. Sk3 = Run.Sk.nk + f - D * nk
    @. nk = nk0 + dt * Sk3


    St_k!(Run)
    @. Sk4 = Run.Sk.nk + f - D * nk

    @. nk = nk0 + dt * (Sk1 + 2.0 * Sk2 + 2.0 * Sk3 + Sk4) / 6.0
    # @. nk = ( nk -.5 * dt * D * nk)/(1. + .5 * dt * D)
    clean_waveaction!(Run.Nk)

    Run.t += dt
    return nothing
end


function init_temporal_scheeme!(temp_scheeme::ETD4_step, Run, dt)
    println("Using ETD4 temporal scheeme")
    println("")

    lk1 = temp_scheeme.arrays.lk1
    lk2 = temp_scheeme.arrays.lk2
    cnl1 = temp_scheeme.arrays.cnl1
    cnl2 = temp_scheeme.arrays.cnl2
    cnl3 = temp_scheeme.arrays.cnl3
    cnl4 = temp_scheeme.arrays.cnl4
    fs1 = temp_scheeme.arrays.fs1

    D = Run.FD.D


    cc = -D
    ch = dt * cc
    ch_2 = ch / 2.0
    c3h2 = cc .^ 3 * dt^2
    @. lk1 = exp(ch)
    @. lk2 = exp(ch_2)
    #    @. cnl1 = expm1(ch_2) / cc
    nl1 = 0.5 * dt * map(x -> abs(x) > 5.e-2 ? expm1(x) / x : 1 + x / 2 + x^2 / 6 + x^3 / 24 + x^4 / 120 + x^5 / 720 + x^6 / 5040 + x^7 / 40320 + x^8 / 362880, ch_2)
    nl2 = dt * map(x -> abs(x) > 5.e-3 ? (-4.0 - x + exp(x) * (4.0 - 3.0 * x + x^2)) / x^3 : 1.0 / 6.0 + x / 6.0 + 3.0 * x^2 / 40.0 + x^3 / 45.0 + 5.0 * x^4 / 1008.0, ch)
    nl3 = dt * map(x -> abs(x) > 5.e-3 ? 2.0 * (2.0 + x + exp(x) * (-2.0 + x)) / x^3 : 2.0 * (1.0 / 6.0 + x / 12.0 + x^2 / 40.0 + x^3 / 180.0 + x^4 / 1008.0), ch)
    nl4 = dt * map(x -> abs(x) > 5.e-3 ? (-4.0 - 3.0 * x - x^2 + exp(x) * (4.0 - x)) / x^3 : 1.0 / 6.0 - x^2 / 120.0 - x^3 / 360.0 - x^4 / 1680.0, ch)
    @. cnl1 = nl1
    @. cnl2 = nl2
    @. cnl3 = nl3
    @. cnl4 = nl4
    return nothing
end

@doc raw"""
    advance!(::ETD4_step, Run, dt)

Make one time step with `ETD4_step`.
* `Run`: run structure
* `dt`: time step

"""
function advance!(temp_scheeme::ETD4_step, Run, dt)

    nk = Run.Nk.nk
    nk0 = temp_scheeme.arrays.nk0
    f = Run.FD.f

    lk1 = temp_scheeme.arrays.lk1
    lk2 = temp_scheeme.arrays.lk2
    cnl1 = temp_scheeme.arrays.cnl1
    cnl2 = temp_scheeme.arrays.cnl2
    cnl3 = temp_scheeme.arrays.cnl3
    cnl4 = temp_scheeme.arrays.cnl4

    Sk1 = temp_scheeme.arrays.Sk1
    Sk2 = temp_scheeme.arrays.Sk2
    Sk3 = temp_scheeme.arrays.Sk3
    Sk4 = temp_scheeme.arrays.Sk4
    fs1 = temp_scheeme.arrays.fs1

    @. nk0 = Run.Nk.nk

    St_k!(Run)
    @. Sk1 = Run.Sk.nk
    @. nk = lk2 * nk0 + cnl1 * (Sk1 + f)
    @. fs1 = nk
    #    @threads for n in eachindex(nk)
    #        nk[n] =lk[n] * nk[n] +cnl1[n]* (Sk[n] +  f[n]) 
    #    end
    clean_waveaction!(Run.Nk)

    St_k!(Run)
    @. Sk2 = Run.Sk.nk
    @. nk = lk2 * nk0 + cnl1 * (Sk2 + f)
    #    @threads for n in eachindex(nk)
    #        nk[n] = nk[n] +cnl2[n]* (Sk[n] - Sktmp[n]) 
    #    end
    clean_waveaction!(Run.Nk)

    St_k!(Run)
    @. Sk3 = Run.Sk.nk
    @. nk = lk2 * fs1 + cnl1 * (2.0 * Sk3 - Sk1 + f)
    #    @threads for n in eachindex(nk)
    #        nk[n] = nk[n] +cnl2[n]* (Sk[n] - Sktmp[n]) 
    #    end
    clean_waveaction!(Run.Nk)

    St_k!(Run)
    @. Sk4 = Run.Sk.nk

    @. nk = lk1 * nk0 + cnl2 * (Sk1 + f) + cnl3 * (Sk2 + Sk3 + 2.0 * f) + cnl4 * (Sk4 + f)
    #    @threads for n in eachindex(nk)
    #        nk[n] = nk[n] +cnl2[n]* (Sk[n] - Sktmp[n]) 
    #    end
    clean_waveaction!(Run.Nk)

    Run.t += dt
    return nothing
end



@doc raw"""
    get_T_nonlinear(Run)

Non linear time based on the wave action spectrum

``\tau_{\rm nl} = \frac{1}{\max\limits_{\bf k} \left| \frac{St_{\bf k}}{n_{\bf k}} \right| }.``

We use only the energetic modes with ``n_{\bf k} > 10^{-50}``.

Note: We use ``\tau_{\rm nl}`` to fix the adaptative time step (see [`adaptative_time_step`](@ref)). 
"""
function get_T_nonlinear(Run)
    inds = Run.Nk.nk .> 1.e-50
    if any(inds)
        @. Run.F1.nk[inds] = Run.Sk.nk[inds] ./ Run.Nk.nk[inds]
        return 1.0 / maximum(abs.(Run.F1.nk[inds]))
    else
        return -1.0
    end
    # @. Run.F1.nk = Run.Sk.nk ./ (1.e-40 + Run.Nk.nk)

    # return 1.0 / maximum(Run.F1.nk)
end

@doc raw"""
    adaptative_time_step(Run, dtmin, dtmax, dt, cmin=0.05, cmax=0.5)

Return the adaptated time step based on the wave action spectrum.  
* `Run`: run structure.   
* `dtmin`: minimal time step.   
* `dtmax`: maximal time step.   
* `dt`: current time step.   
* `cmin`: lower threshold for the ratio of the nonlinear time step and current time step.    
* `cmax`: upper threshold for the ratio of the nonlinear time step and current time step.  
See also [`get_T_nonlinear`](@ref).

!!! warning 
    This adaptative time step may not work for wave action that decrease very quickly.
"""
function adaptative_time_step(Run, dtmin, dtmax, dt, cmin=0.05, cmax=0.5)
    τnl = get_T_nonlinear(Run)
    if τnl > 0
        if ((dt / τnl) > cmax && dt > dtmin) || dt > dtmax# time step too large
            dt = dt / 1.25
            println(" Time step decreased. dt = ", dt)
        end
        if ((dt / τnl) < cmin && dt < dtmax) || dt < dtmin # time step too large
            dt = 1.25 * dt
            println(" Time step increased. dt = ", dt)
        end
    end
    return dt
end