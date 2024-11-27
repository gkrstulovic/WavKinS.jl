push!(LOAD_PATH, "../")
using WavKinS
using BenchmarkTools

## Defining the simulation for Acoustic2D
n = 8
M = 2^n
kmin = 1e-4
kmax = 1e+0
Nk = wave_spectrum(kmin,kmax,M)

# The Acoustic2D structure contains defines all the methods for integration, see help ?Acoustic2D
Run = Acoustic2D(Nk; interp_scheeme=WavKinS.lin_interp, time_stepping_scheeme=WavKinS.Euler_step)

println(" ********************************************************************** ")
println(" Collisional integral benckmark with M=",M," using  ",Threads.nthreads(), " threads")
@btime WavKinS.St_k!(Run)
@btime advance!(Run.time_stepping, Run, 0.)
println(" ********************************************************************** ")

## Defining the simulation for Stratified_Asymp

M = 80
Mh = M
Mz = M
khmin = 1e-3
khmax = 1e0
kzmin = 1e-3
kzmax = 1e0
Nk = wave_spectrum_khkz(khmin, khmax, Mh, kzmin, kzmax, Mz);

# The Stratified_Asymp structure contains defines all the methods for integration, see help ?Stratified_Asymp
Run = Stratified_Asymp(Nk; interp_scheeme=WavKinS.bilin_interp_khkz, time_stepping_scheeme=WavKinS.RK2_step);

println(" ********************************************************************** ")
println(" Collisional integral benckmark with M=",M," using  ",Threads.nthreads(), " threads")
@btime WavKinS.St_k!(Run)
@btime advance!(Run.time_stepping, Run, 0.)
println(" ********************************************************************** ")
