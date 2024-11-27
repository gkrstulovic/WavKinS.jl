push!(LOAD_PATH, "../")
using WavKinS
using BenchmarkTools
 using Profile
 using PProf

PowerLawField = false
kmin = 1e-2
kmax = 1e1
M = 1600


if PowerLawField
    Nk = waveaction(kmin,kmax,M)
    Run = Acoustic2D(Nk)
else
    Nk = wave_spectrum(kmin,kmax,M)
    Run = Acoustic2D(Nk;interp_scheeme=WavKinS.BS_interp)
end

#@btime WavKinS.St_k!(Run)

#@btime WavKinS.St_k_threaded!(Run)

##
#WavKinS.St_k!(Run)
#@profile WavKinS.St_k!(Run)

# WavKinS.St_k!(Run)
# WavKinS.St_k!(Run)

WavKinS.St_k!(Run)

Profile.Allocs.clear()
Profile.Allocs.@profile sample_rate=0.001 WavKinS.St_k!(Run)
results = Profile.Allocs.fetch()

PProf.Allocs.pprof(results)

# @profile WavKinS.St_k_threaded!(Run)
# Export pprof profile and open interactive profiling web interface.
# pprof()