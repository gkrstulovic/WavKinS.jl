using Base.Threads
# Here we define collisional integrals for Acoustic runs


function St_k!(Run::Bogoliubov3D)
    Nk = Run.Nk
    nk = Nk.nk
    Sk = Run.Sk.nk

    ξ=Run.ξ
    c=Run.c
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
                    k2 = k2_Bogoliubov(k,k1;ξ=ξ)
                    nk2 = val_nk(interp_scheeme, Nk, k2)
                    S123 = V123squared_Bogo3D(k, k1, k2; V0=V0GP,ξ=ξ) #Vsquared_Bogoliubov_3D
                    S123 *= θkmax(k1,kmax)*θkmax(k2,kmax)
                    F1[j] = S123 * nk[j] * nk2 - S123 * (nk[j] + nk2) * nk[i]
                    F1[j] = F1[j]/Det_Bogo(k2;c=c,ξ=ξ)
                end
                Sk[i] += integrate_with_log_bins(Run.FSt[id], 1, i)

                F1 .= 0.0
                for j = i:M
                    k1 = kk[j]
                    k2 = k2_Bogoliubov(k,k1;ξ=ξ)
                    nk2 = val_nk(interp_scheeme, Nk, k2)
                    S123 = V123squared_Bogo3D(k1, k, k2; V0=V0GP,ξ=ξ)
                    S123 *= θkmax(k1,kmax)*θkmax(k2,kmax)
                    F1[j] = -2 * S123 * nk2 * nk[i] + 2 * S123 * nk[j] * nk[i] + 2 * S123 * nk[j] * nk2
                    F1[j] = F1[j]/Det_Bogo(k2;c=c,ξ=ξ)
                end
                Sk[i] += integrate_with_log_bins(Run.FSt[id], i + 1, M)

                Sk[i] *= Sfactor / k^2

            end
        end
    end
end
