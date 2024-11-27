using Base.Threads
# Here we define collisional integrals for Smoluchowski run

function St_k!(Run::Smoluchowski)
    Nk = Run.Nk
    nk = Nk.nk
    Sk = Run.Sk.nk
    K = Run.K

    kk = Nk.kk
    M = Nk.M

    interp_scheeme = Run.interp_scheeme
    update_coeff_interp!(interp_scheeme, Nk)

    @sync for id = 1:Threads.nthreads()
        inds = Run.partition[id]
        F1 = Run.FSt[id].nk
        @spawn begin
            for i in inds
                k = kk[i]
                Sk[i] = 0.0

                # First (source) integral
                F1 .= 0.0
                for j = 1:i
                    p = kk[j]
                    q = k - p
                    np = val_nk(interp_scheeme, Nk, p)
                    nq = val_nk(interp_scheeme, Nk, q)
                    F1[j] = K(p,q) * np * nq 
                end
                Sk[i] += integrate_with_log_bins(Run.FSt[id], 1, i) / 2

                # Second (sink) integral
                F1 .= 0.0
                for j = 1:M
                    p = kk[j]
                    np = val_nk(interp_scheeme, Nk, p)
                    F1[j] = - K(k,p) * nk[i] * np
                end
                Sk[i] += integrate_with_log_bins(Run.FSt[id], 1, M)

            end
        end
    end
end
