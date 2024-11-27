using Base.Threads
# Here we define collisional integrals for Acoustic runs

function St_k!(Run::Acoustic2D)
    Nk = Run.Nk
    nk = Nk.nk
    Sk = Run.Sk.nk

    kk = Nk.kk
    M = Nk.M

    interp_scheeme = Run.interp_scheeme
    update_coeff_interp!(interp_scheeme, Nk)

    Sfactor = 2 * pi / (sqrt(6) * Run.a * Run.c)
    V0GP = 3 * sqrt(Run.c / 2.0) / 4.0

    @sync for id = 1:Threads.nthreads()
        inds = Run.partition[id]
        F1 = Run.FSt[id].nk
        @spawn begin
            for i in inds
                k = kk[i]
                Sk[i] = 0.0

                F1 .= 0.0
                for j = 1:i
                    k1 = kk[j]
                    k2 = k - k1
                    nk2 = val_nk(interp_scheeme, Nk, k2)
                    S123 = S_2DAcoustic(k, k1, k2; V0=V0GP)
                    F1[j] = S123 * nk[j] * nk2 - S123 * (nk[j] + nk2) * nk[i]
                end
                #  Sk[i] += integrate(Run.FSt[id], 1, i)
                Sk[i] += integrate_with_log_bins(Run.FSt[id], 1, i)


                F1 .= 0.0
                for j = i:M
                    k1 = kk[j]
                    k2 = k1 - k
                    nk2 = val_nk(interp_scheeme, Nk, k2)
                    S123 = S_2DAcoustic(k, k1, k2; V0=V0GP)
                    F1[j] = -2 * S123 * nk2 * nk[i] + 2 * S123 * nk[j] * nk[i] + 2 * S123 * nk[j] * nk2
                end
                # Sk[i] += integrate(Run.FSt[id], i + 1, M)
                Sk[i] += integrate_with_log_bins(Run.FSt[id], i + 1, M)

                Sk[i] *= Sfactor / k

            end
        end
    end
end

function St_k!(Run::Acoustic3D)
    Nk = Run.Nk
    nk = Nk.nk
    Sk = Run.Sk.nk

    kk = Nk.kk
    M = Nk.M

    interp_scheeme = Run.interp_scheeme
    update_coeff_interp!(interp_scheeme, Nk)

    Sfactor = 4 * pi^2 / (Run.c)
    V0GP = 3 * sqrt(Run.c / 2.0) / 4.0
    kmax = kk[end]

    @sync for id = 1:Threads.nthreads()
        inds = Run.partition[id]
        F1 = Run.FSt[id].nk
        @spawn begin
            for i in inds
                k = kk[i]
                Sk[i] = 0.0

                F1 .= 0.0
                for j = 1:i
                    k1 = kk[j]
                    k2 = k - k1
                    nk2 = val_nk(interp_scheeme, Nk, k2)
                    S123 = S_3DAcoustic(k, k1, k2; V0=V0GP)
                    S123 *= θkmax(k1,kmax)*θkmax(k2,kmax)
                    F1[j] = S123 * nk[j] * nk2 - S123 * (nk[j] + nk2) * nk[i]
                end
                #  Sk[i] += integrate(Run.FSt[id], 1, i)
                Sk[i] += integrate_with_log_bins(Run.FSt[id], 1, i)

                F1 .= 0.0
                for j = i:M
                    k1 = kk[j]
                    k2 = k1 - k
                    nk2 = val_nk(interp_scheeme, Nk, k2)
                    S123 = S_3DAcoustic(k, k1, k2; V0=V0GP)
                    S123 *= θkmax(k1,kmax)*θkmax(k2,kmax)
                    F1[j] = -2 * S123 * nk2 * nk[i] + 2 * S123 * nk[j] * nk[i] + 2 * S123 * nk[j] * nk2
                end
                # Sk[i] += integrate(Run.FSt[id], i + 1, M)
                Sk[i] += integrate_with_log_bins(Run.FSt[id], i + 1, M)

                Sk[i] *= Sfactor / k^2

            end
        end
    end
end
