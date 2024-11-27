using Base.Threads
# Here we define collisional integrals for MMT runs

function St_k!(Run::NLS3D; compute_Sk=true, compute_γk=false, compute_ηk=false)
    Nk = Run.Nk
    nk = Nk.nk
    Sk = Run.Sk.nk
    γk = Run.γk
    ηk = Run.ηk

    logλ = Run.Nk.logλ
    logλω = Run.FSt2[1].logλ
    ωω = Run.FSt1[1].kk

    kk = Nk.kk
    M = Nk.M
    kmax = kk[end]
    ωmax = kmax^2

    interp_scheeme = Run.interp_scheeme
    update_coeff_interp!(interp_scheeme, Nk)

    Sfactor = 4 * pi^3

    @sync for id = 1:Threads.nthreads()
        inds = Run.partition[id]
        dF3 = Run.FSt[id].nk
        dF2 = Run.FSt1[id].nk
        dγk3 = Run.FSt2[id].nk
        dηk3 = Run.FSt3[id].nk
        dγk2 = Run.FSt4[id].nk
        dηk2 = Run.FSt5[id].nk
        @spawn begin
            for i in inds
                k = kk[i]
                ω = ωω[i]
                Sk[i] = 0.0

                dF2 .= 0.0
                dγk2 .= 0.0
                dηk2 .= 0.0
                for j2 = 1:i # integration is perfomed in the ω2-ω3 plane (standard convention)
                    ω2 = ωω[j2]
                    dF3 .= 0.0
                    dγk3 .= 0.0
                    dηk3 .= 0.0

                    ω3inf = ω - ω2
                    if j2 == i
                        iω3min = 1
                    else
                        iω3min = ceil(Int, 1 + log(ω3inf / ωω[1]) / logλω) # Mesh point above ω3inf
                        iω3min = max(min(iω3min, M), 1)
                    end
                    for j3 = iω3min:M
                        ω3 = ωω[j3]
                        ω1 = ω2 + ω3 - ω
                        n1 = val_nk(interp_scheeme, Nk, sqrt(ω1))
                        n2 = nk[j2]
                        n3 = nk[j3]
                        T1234squared = T1234squared_NLS3D_in_ω(ω, ω1, ω2, ω3)
                        T1234squared *= θkmax(ω1, ωmax) * θkmax(ω2, ωmax) * θkmax(ω3, ωmax)
                        Det = 1.0

                        if compute_Sk
                            if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                                dγk3[j3] = -T1234squared * (n2 * n3 - n1 * n3 - n1 * n2) / Det
                                dηk3[j3] = T1234squared * (n1 * n2 * n3) / Det
                                @assert dηk3[j3] >= 0 "I: neg vals dηk3 = $dηk3"
                            else
                                dF3[j3] = T1234squared * (n1 * n2 * n3 + nk[i] * n2 * n3 - nk[i] * n1 * n3 - nk[i] * n1 * n2) / Det
                            end
                        else
                            if compute_γk
                                dγk3[j3] = -T1234squared * (n2 * n3 - n1 * n3 - n1 * n2) / Det
                            end
                            if compute_ηk
                                dηk3[j3] = T1234squared * (n1 * n2 * n3) / Det
                            end
                        end

                    end
                    if compute_Sk
                        if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                            dγk2[j2] = integrate_with_log_bins(Run.FSt2[id], iω3min + 1, M)
                            dγk2[j2] += integrate_with_log_bins_segment(0.0, dγk3[iω3min], ω3inf, ωω[iω3min])# add missing part of the integral           

                            dηk2[j2] = integrate_with_log_bins(Run.FSt3[id], iω3min + 1, M)
                            dηk2[j2] += integrate_with_log_bins_segment(0.0, dηk3[iω3min], ω3inf, ωω[iω3min])# add missing part of the integral           

                        else
                            dF2[j2] = integrate_with_log_bins(Run.FSt[id], iω3min + 1, M)
                            dF2[j2] += integrate_with_log_bins_segment(0.0, dF3[iω3min], ω3inf, ωω[iω3min])# add missing part of the integral           
                        end
                    else
                        if compute_γk
                            dγk2[j2] = integrate_with_log_bins(Run.FSt2[id], iω3min + 1, M)
                            dγk2[j2] += integrate_with_log_bins_segment(0.0, dγk3[iω3min], ω3inf, ωω[iω3min])# add missing part of the integral           
                        end
                        if compute_ηk
                            dηk2[j2] = integrate_with_log_bins(Run.FSt3[id], iω3min + 1, M)
                            dηk2[j2] += integrate_with_log_bins_segment(0.0, dηk3[iω3min], ω3inf, ωω[iω3min])# add missing part of the integral           
                        end
                    end
                    # second part

                end

                for j2 = (i+1):M # integration is perfomed in the ω2-ω3 plane (standard convention)
                    ω2 = ωω[j2]
                    dF3 .= 0.0
                    dγk3 .= 0.0
                    dηk3 .= 0.0


                    ω3max = ω + ωmax - ω2
                    iω3max = ceil(Int, 1 + log(ω3max / ωω[1]) / logλω) - 1# Mesh point below ω3max

                    for j3 = 1:iω3max
                        ω3 = ωω[j3]
                        ω1 = ω2 + ω3 - ω
                        n1 = val_nk(interp_scheeme, Nk, sqrt(ω1))
                        n2 = nk[j2]
                        n3 = nk[j3]
                        T1234squared = T1234squared_NLS3D_in_ω(ω, ω1, ω2, ω3)
                        T1234squared *= θkmax(ω1, ωmax) * θkmax(ω2, ωmax) * θkmax(ω3, ωmax)
                        Det = 1.0

                        if compute_Sk
                            if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                                dγk3[j3] = -T1234squared * (n2 * n3 - n1 * n3 - n1 * n2) / Det
                                dηk3[j3] = T1234squared * (n1 * n2 * n3) / Det
                                @assert dηk3[j3] >= 0 "I: neg vals dηk = $dηk"
                            else
                                dF3[j3] = T1234squared * (n1 * n2 * n3 + nk[i] * n2 * n3 - nk[i] * n1 * n3 - nk[i] * n1 * n2) / Det
                            end
                        else
                            if compute_γk
                                dγk3[j3] = -T1234squared * (n2 * n3 - n1 * n3 - n1 * n2) / Det
                            end
                            if compute_ηk
                                dηk3[j3] = T1234squared * (n1 * n2 * n3) / Det
                            end
                        end

                    end

                    if iω3max > 0 # 
                        if compute_Sk
                            if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                                dγk2[j2] = integrate_with_log_bins(Run.FSt2[id], 1, iω3max)
                                dγk2[j2] += integrate_with_log_bins_segment(dγk3[iω3max], 0.0, ωω[iω3max], ω3max)# add missing part of the integral   

                                dηk2[j2] = integrate_with_log_bins(Run.FSt3[id], 1, iω3max)
                                dηk2[j2] += integrate_with_log_bins_segment(dηk3[iω3max], 0.0, ωω[iω3max], ω3max)# add missing part of the integral   

                            else
                                dF2[j2] = integrate_with_log_bins(Run.FSt[id], 1, iω3max)
                                dF2[j2] += integrate_with_log_bins_segment(dF3[iω3max], 0.0, ωω[iω3max], ω3max)# add missing part of the integral   
                            end
                        else
                            if compute_γk
                                dγk2[j2] = integrate_with_log_bins(Run.FSt2[id], 1, iω3max)
                                dγk2[j2] += integrate_with_log_bins_segment(dγk3[iω3max], 0.0, ωω[iω3max], ω3max)# add missing part of the integral   
                            end
                            if compute_ηk
                                dηk2[j2] = integrate_with_log_bins(Run.FSt3[id], 1, iω3max)
                                dηk2[j2] += integrate_with_log_bins_segment(dηk3[iω3max], 0.0, ωω[iω3max], ω3max)# add missing part of the integral   
                            end
                        end
                    end

                end

                if compute_Sk
                    if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                        γk[i] = integrate_with_log_bins( Run.FSt4[id])
                        γk[i] *= Sfactor / k
                        ηk[i] = integrate_with_log_bins(Run.FSt5[id])
                        ηk[i] *= Sfactor / k
                        Sk[i] = -γk[i] * nk[i] + ηk[i]
                    else
                        Sk[i] = integrate_with_log_bins(Run.FSt1[id])
                        Sk[i] *= Sfactor / k
                    end
                else
                    if compute_γk
                        γk[i] =integrate_with_log_bins( Run.FSt4[id])
                        γk[i] *= Sfactor / k
                    end
                    if compute_ηk
                        ηk[i] = integrate_with_log_bins( Run.FSt5[id])
                        ηk[i] *= Sfactor / k
                    end
                end



            end
        end
    end
end