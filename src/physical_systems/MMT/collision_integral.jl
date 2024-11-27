using Base.Threads
# Here we define collisional integrals for MMT runs

function St_k!(Run::MMT; compute_Sk=true, compute_γk=false, compute_ηk=false)
    Nk = Run.Nk
    nk = Nk.nk
    Sk = Run.Sk.nk
    γk = Run.γk
    ηk = Run.ηk


    kk = Nk.kk
    M = Nk.M
    kmax = kk[end]

    interp_scheeme = Run.interp_scheeme
    update_coeff_interp!(interp_scheeme, Nk)

    Sfactor = 4pi

    @sync for id = 1:Threads.nthreads()
        inds = Run.partition[id]
        F1 = Run.FSt[id].nk
        dγk = Run.FSt1[id].nk
        dηk = Run.FSt2[id].nk
        @spawn begin
            for i in inds
                k = kk[i]
                Sk[i] = 0.0

                F1 .= 0.0
                dγk .= 0.0
                dηk .= 0.0
                for j = 1:M # we integrate from -∞ to 0, so k3 should be negative
                    if j != i
                        k3 = -kk[j] #
                        k1 = get_k1_MMT(k, k3)
                        k2 = k + k1 - k3
                        n3 = nk[j]
                        n2 = val_nk(interp_scheeme, Nk, abs(k2))
                        n1 = val_nk(interp_scheeme, Nk, abs(k1))
                        T1234squared = T1234squared_MMT(k, k1, k2, k3; β=Run.β)
                        T1234squared *= θkmax(k1, kmax) * θkmax(k2, kmax) * θkmax(k3, kmax)
                        Det = Det_MMT(k, k1, k2, k3)

                        if compute_Sk
                            if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                                dγk[j] = T1234squared * (n1 * n3 + n1 * n2 - n2 * n3) / Det
                                dηk[j] = T1234squared * (n1 * n2 * n3) / Det
                                @assert dηk[j] >= 0 "I: neg vals dηk = $dηk"
                            else
                                F1[j] = T1234squared * (n1 * n2 * n3 + nk[i] * n2 * n3 - nk[i] * n1 * n3 - nk[i] * n1 * n2) / Det
                            end
                        else
                            if compute_γk
                                dγk[j] = T1234squared * (n1 * n3 + n1 * n2 - n2 * n3) / Det
                            end
                            if compute_ηk
                                dηk[j] = T1234squared * (n1 * n2 * n3) / Det
                            end
                        end
                    end
                end
                if compute_Sk
                    if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                        γk[i] = integrate_with_log_bins(Run.FSt1[id], 1, M)
                        ηk[i] = integrate_with_log_bins(Run.FSt2[id], 1, M)
                    else
                        Sk[i] = integrate_with_log_bins(Run.FSt[id], 1, M)
                    end
                else
                    if compute_γk
                        γk[i] = integrate_with_log_bins(Run.FSt1[id], 1, M)
                    end
                    if compute_ηk
                        ηk[i] = integrate_with_log_bins(Run.FSt2[id], 1, M)
                    end
                end

                F1 .= 0.0
                dγk .= 0.0
                dηk .= 0.0
                for j = 1:M # we integrate from 0 to ∞
                    if j != i
                        k3 = kk[j] #
                        k1 = get_k1_MMT(k, k3)
                        k2 = k + k1 - k3
                        n3 = nk[j]
                        n2 = val_nk(interp_scheeme, Nk, abs(k2))
                        n1 = val_nk(interp_scheeme, Nk, abs(k1))
                        T1234squared = T1234squared_MMT(k, k1, k2, k3; β=Run.β)
                        T1234squared *= θkmax(k1, kmax) * θkmax(k2, kmax) * θkmax(k3, kmax)
                        Det = Det_MMT(k, k1, k2, k3)

                        if compute_Sk
                            if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                                dγk[j] = T1234squared * (n1 * n3 + n1 * n2 - n2 * n3) / Det
                                dηk[j] = T1234squared * (n1 * n2 * n3) / Det
                                @assert dηk[j] >= 0 "I: neg vals dηk = $dηk"
                            else
                                F1[j] = T1234squared * (n1 * n2 * n3 + nk[i] * n2 * n3 - nk[i] * n1 * n3 - nk[i] * n1 * n2) / Det
                            end
                        else
                            if compute_γk
                                γk[j] = T1234squared * (n1 * n3 + n1 * n2 - n2 * n3) / Det
                            end
                            if compute_ηk
                                dηk[j] = T1234squared * (n1 * n2 * n3) / Det
                            end
                        end
                    end
                end

                if compute_Sk
                    if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                        γk[i] += integrate_with_log_bins(Run.FSt1[id], 1, M)
                        γk[i] *= Sfactor
                        ηk[i] += integrate_with_log_bins(Run.FSt2[id], 1, M)
                        ηk[i] *= Sfactor
                        Sk[i] = -γk[i] * nk[i] + ηk[i]
                    else
                        Sk[i] += integrate_with_log_bins(Run.FSt[id], 1, M)
                        Sk[i] *= Sfactor
                    end
                else
                    if compute_γk
                        γk[i] += integrate_with_log_bins(Run.FSt1[id], 1, M)
                        γk[i] *= Sfactor
                    end
                    if compute_ηk
                        ηk[i] += integrate_with_log_bins(Run.FSt2[id], 1, M)
                        ηk[i] *= Sfactor
                    end
                end



            end
        end
    end
end
