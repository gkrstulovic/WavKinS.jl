using WavKinS

## Defining the wave_spectrum structure

kmin = 1e-3;
kmax = 1e0;

for M ∈ 2 .^ (5:12)
    Nk = wave_spectrum(kmin, kmax, M)

    ## Creating the Run Acoustic2D, filling the waveaction and computing the collisional integral
    Run = Acoustic2D(Nk)
    kk = Run.Nk.kk # we rename the wavevector mesh for shorter use
    nk = Run.Nk.nk # we rename the wavefunction mesh for shorter use
    @. nk = 1.0 / Run.ω(kk) # the @. is the julia macro for pointwise operations


    ## Compute the collisional integral and print th
    St_k!(Run) # compute the collision integral

    ## Compute the L2-norm of the collisional integral

    @. Run.Sk.nk = 2pi * kk * Run.Sk.nk^2 # we use the same waveaction structure to compute the square
    L2norm = sqrt(integrate_with_log_bins(Run.Sk)) # we integrate collisional integral squared using WavKinS integration

    println("M =", M, " L2-norm = ", L2norm)

end