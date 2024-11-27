# Tests for time_stepping

push!(LOAD_PATH, "../")
using WavKinS
using GLMakie
using TimerOutputs


function nk_theo(t,d,f,n0)
    c = sqrt(Complex(4*f - d^2))
    nk = ( d + c * tan( atan((-d+2*n0)/c) + c * t / 2) ) / 2
    return real(nk)
end



## Test collision_integrals

println("---------------------------------------------------------------------")
println("Error of time stepping method")
println("")



## Defining the simulation
Tfinal = 1.0; # Final time of the simulation
dtplot = 1.0; # Plot every tplot times

M = 64
kmin = 1e-2
kmax = 1e0
Nk = wave_spectrum(kmin,kmax,M)
Run = TestTimeStepping(Nk; time_stepping_scheeme=WavKinS.RK4_step) # Select time stepping method here

## Setting forcing, dissipation and initial condition            
kk = Nk.kk
kd = kmax / 2

n0 = kk;
@. Run.FD.f = 1.0
@. Run.FD.D = (kk / kd)^2


dts = Vector{Float64}()
errors = Vector{Float64}()

for M = 0:10

    dt = 0.1 / 2^M; # time step of the simulation

    Run.t = 0.0
    @. Run.Nk.nk = n0;

    init_temporal_scheeme!(Run.time_stepping, Run, dt)
    
    @time while Run.t <= Tfinal

        advance!(Run.time_stepping, Run, dt)

        if mod(Run.t,dtplot) < dt
            fig, ax = plot_wave_action!(Run)
            plot_theo!(ax, kk, nk_theo.(Run.t,Run.FD.D,Run.FD.f,n0))
            display(fig)
            sleep(1)
        end
    end

    nk = Run.Nk.nk
    nk_th = nk_theo.(Run.t,Run.FD.D,Run.FD.f,n0)
    tmp = sum( sqrt.(  (nk - nk_th).^2 )  )
    append!(dts,dt)
    append!(errors,tmp)
end


fig, ax = scatter(dts,errors,label="Time_Stepping")
ylims!(1e-10, 1e1)
ax.xlabel = "dt"
ax.xscale = log10
ax.ylabel = "Error"
ax.yscale = log10
ax.title = "Error of time stepping"
lines!(ax, dts, 1e2 * dts.^(1); label=L"\mathrm{d}t^{1}", color="red", linestyle=:dash)
lines!(ax, dts, 1e2 * dts.^(2); label=L"\mathrm{d}t^{2}", color="magenta", linestyle=:dash)
lines!(ax, dts, 1e2 * dts.^(4); label=L"\mathrm{d}t^{4}", color="blue", linestyle=:dash)
axislegend(ax, position = :lt)
display(fig)
sleep(5)

println("")
println("---------------------------------------------------------------------")


