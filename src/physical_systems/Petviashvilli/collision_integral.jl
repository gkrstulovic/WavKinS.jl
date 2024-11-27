using Base.Threads
# Here we define collisional integrals for Petviashvilli runs


function St_k!(Run::Petviashvilli)
    Nk = Run.Nk
    nk = Nk.nk
    Sk = Run.Sk.nk

    kkx = Nk.kkh
    kky = Nk.kkz
    Mx = Nk.Mh
    My = Nk.Mz

    interp_scheeme = Run.interp_scheeme
    update_coeff_interp!(interp_scheeme, Nk)

    Sfactor = 2 * pi

    @sync for id = 1:Threads.nthreads()
        inds = Run.partition[id]
        F1 = Run.FSt[id].nk
        Mk1 = Run.FSt[id].M
        kk1 = Run.FSt[id].kk
        logλ1 = Run.FSt[id].logλ
        @spawn begin
            #@inbounds for iy in eachindex(kky), ix in eachindex(kkx)
            @inbounds for iy in inds[2], ix in inds[1]

                kx = kkx[ix]
                ky = kky[iy]
                Sk[ix, iy] = 0.0

                # First integral
                F1 .= 0.0
                #i1x_kx = max(1, 1 - Ik1[1] + ceil(Int, log(kx) / logλ1))
                i1x_kx = ceil(Int, 1+log(kx/kk1[1]) / logλ1) # TODO: Remove previous line if ok
                i1x_kx = min(i1x_kx, Mk1)
                i1x_kx = ix # TODO: What is the purpose of the 3 previous lines?
                @inbounds for i1x = 1:(i1x_kx-1) # no contibution for i1x=ix as k2x=0
                    k1x = kk1[i1x]
                    k2x = kx - k1x
                    k1y_plus, k1y_minus = get_k1y_Petviashvilli(kx, ky, k1x)
                    k2y_plus = abs(ky - k1y_plus)
                    k2y_minus = abs(ky - k1y_minus) # symmetry in ky-> -ky assumed
                    Det = Det_Petviashvilli(kx, ky, k1x)
                    n1_plus = val_nk(interp_scheeme, Nk, k1x, abs(k1y_plus))
                    n1_minus = val_nk(interp_scheeme, Nk, k1x, abs(k1y_minus))
                    n2_plus = val_nk(interp_scheeme, Nk, k2x, abs(k2y_plus))
                    n2_minus = val_nk(interp_scheeme, Nk, k2x, abs(k2y_minus))
                    S123 = S_Petviashvilli(kx, k1x, k2x)
                    F1[i1x] = S123 * (n1_plus * n2_plus - n1_plus * nk[ix, iy] - n2_plus * nk[ix, iy]) / Det
                    F1[i1x] += S123 * (n1_minus * n2_minus - n1_minus * nk[ix, iy] - n2_minus * nk[ix, iy]) / Det
                end
                Sk[ix, iy] = integrate_with_log_bins(Run.FSt[id], 1, i1x_kx)
                #Sk[ix, iy] += integrate_with_log_bins_segment(F1[i1x_kx-1], 0, kk1[i1x_kx-1], kx) # add missing segment between last grid point and kx

                #Second integral
                F1 .= 0.0
                if kx <= ky / sqrt(3) # simple part, no singularities
                    #   println("In simple part")
                    ik1end = Mk1
                else
                    #  println("In singular part")
                    kxsup = get_ksup_Petviashvilli(kx, ky)
                    #ikMax = max(1, 1 - Ik1[1] + ceil(Int, max(log(kxsup) / logλ1))) # finding index such that λ1^ikMax∼ k
                    ikMax = max(1,ceil(Int, 1+log(kxsup/kk1[1]) / logλ1)) # TODO: Remove previous line if ok
                    ik1end = min(ikMax - 1, Mk1)

                end

                F1tmp = 0.0
                F1ini = 0.0
                @inbounds for i1x = (i1x_kx+1):ik1end
                    k1x = kk1[i1x]
                    k2x = k1x - kx
                    k1y_plus, k1y_minus = get_k1y_Petviashvilli(kx, ky, k1x)
                    k2y_plus = abs(k1y_plus - ky)
                    k2y_minus = abs(k1y_minus - ky) # symmetry in ky-> -ky assumed
                    Det = Det_Petviashvilli(kx, ky, k1x)
                    n1_plus = val_nk(interp_scheeme, Nk, k1x, abs(k1y_plus))
                    n1_minus = val_nk(interp_scheeme, Nk, k1x, abs(k1y_minus))
                    n2_plus = val_nk(interp_scheeme, Nk, k2x, abs(k2y_plus))
                    n2_minus = val_nk(interp_scheeme, Nk, k2x, abs(k2y_minus))
                    S123 = S_Petviashvilli(kx, k1x, k2x)
                    F1tmp = -2 * S123 * (nk[ix, iy] * n2_plus - n1_plus * nk[ix, iy] - n1_plus * n2_plus)
                    F1tmp += -2 * S123 * (nk[ix, iy] * n2_minus - n1_minus * nk[ix, iy] - n1_minus * n2_minus)
                    F1[i1x] = F1tmp / Det
                    # if i1x == (i1x_kx + 1)
                    #     F1ini = F1[i1x]
                    # end
                end
                Sk[ix, iy] += integrate_with_log_bins(Run.FSt[id], i1x_kx + 1, ik1end)
                # if i1x_kx < Mk1
                #     #    Sk[ix, iy] += integrate_with_log_bins_segment(0.0, F1ini, kx, kk1[i1x_kx+1]) # add missing segment between kx and first grid point
                # end

                if kx > ky / sqrt(3) #adding contribution form singulairty
                    Sk[ix, iy] += F1tmp * integrable_singularity_contribution_Petviashvilli(kx, ky, kk1[ik1end], kxsup)
                end
                # Sk[ix, iy] *= Sfactor
            end
        end
    end
    return nothing
end

function St_k!(Run::Petviashvilli_Asymp; compute_Sk=true, compute_γk=false, compute_ηk=false)
    Nk = Run.Nk
    nk = Nk.nk
    Sk = Run.Sk.nk
    γk = Run.γk
    ηk = Run.ηk

    kkx = Nk.kkh
    kky = Nk.kkz
    Mx = Nk.Mh
    My = Nk.Mz

    interp_scheeme = Run.interp_scheeme
    update_coeff_interp!(interp_scheeme, Nk)

    Sfactor = 2 * pi

    @sync for id = 1:Threads.nthreads()
        #    id=1
        inds = Run.partition[id]
        F1 = Run.FSt[id].nk
        dγk = Run.FSt1[id].nk # Attention: it is the same memory than F1
        dηk = Run.FSt2[id].nk
        Mk1 = Run.FSt[id].M
        kk1 = Run.FSt[id].kk
        logλ1 = Run.FSt[id].logλ
        @spawn begin
            #for iy in eachindex(kky), ix in eachindex(kkx)
            @inbounds for iy in inds[2], ix in inds[1]

                kx = kkx[ix]
                ky = kky[iy]
                Sk[ix, iy] = 0.0
                γk[ix, iy] = 0.0
                ηk[ix, iy] = 0.0
                # First integral
                F1 .= 0.0
                dγk .= 0.0
                dηk .= 0.0
                @inbounds for i1y = 1:(iy-1) # In first integral we make k1y-> -k1y to bring the integral to [0 ,ky]. k2x[iy]=0, not computed
                    k1y = kky[i1y]
                    k1x = get_k1x_Petviashvilli_Asymp(kx, ky, -k1y)
                    k2y = ky - (-k1y)
                    k2x = kx - k1x
                    Det = Det_Petviashvilli_Asymp(ky, -k1y)
                    n1 = val_nk(interp_scheeme, Nk, k1x, k1y)
                    n2 = val_nk(interp_scheeme, Nk, k2x, k2y)
                    @assert n1 >= 0 "neg vals n1 = $n1"
                    @assert n2 >= 0 "neg vals n2 = $n2"
                    S123 = S_Petviashvilli(kx, k1x, k2x)
                    if compute_Sk
                        if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                            dγk[i1y] = S123 * (n1 + n2) / Det
                            dηk[i1y] = S123 * (n1 * n2) / Det
                            @assert dηk[i1y] >= 0 "1 neg vals dηk[i1y] = $dηk"
                        else
                            F1[i1y] = S123 * (n1 * n2 - n1 * nk[ix, iy] - n2 * nk[ix, iy]) / Det
                        end
                    else
                        if compute_γk
                            dγk[i1y] = S123 * (n1 + n2) / Det
                        end
                        if compute_ηk
                            dηk[i1y] = S123 * (n1 * n2) / Det
                        end
                    end
                end
                if compute_Sk
                    if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                        γk[ix, iy] = integrate_with_log_bins(Run.FSt1[id], 1, iy)
                        ηk[ix, iy] = integrate_with_log_bins(Run.FSt2[id], 1, iy)
                        if ηk[ix, iy] < 0
                            println("   ")
                            println(" ix=$ix iy=$iy")
                            println(dηk)
                            println(ηk[ix, iy])
                            println("  ")
                        end
                        @assert ηk[ix, iy] >= 0 "1 neg vals ηk[ix, iy] $ix, $iy"
                    else
                        Sk[ix, iy] = integrate_with_log_bins(Run.FSt[id], 1, iy)
                    end
                else
                    if compute_γk
                        γk[ix, iy] = integrate_with_log_bins(Run.FSt1[id], 1, iy)
                    end
                    if compute_ηk
                        ηk[ix, iy] = integrate_with_log_bins(Run.FSt2[id], 1, iy)
                    end
                end

                #i2ky = max(1, 1 - Ik1[1] + ceil(Int, log(2 * ky) / logλ1)) # Mesh point above 2ky
                i2ky = ceil(Int, 1+log(2*ky/kk1[1]) / logλ1) # TODO: Remove previous lines if ok
                i2ky = min(i2ky, My)
                F1[iy] = 0.0
                dγk[iy] = 0.0
                dηk[iy] = 0.0
                @inbounds for i1y = (iy+1):(i2ky-1)
                    k1y = kky[i1y]
                    k1x = get_k1x_Petviashvilli_Asymp(kx, ky, k1y)
                    k2y = ky - k1y
                    k2x = kx - k1x
                    Det = Det_Petviashvilli_Asymp(ky, k1y)

                    n1 = val_nk(interp_scheeme, Nk, k1x, abs(k1y))
                    n2 = val_nk(interp_scheeme, Nk, k2x, abs(k2y))
                    S123 = S_Petviashvilli(kx, k1x, k2x)
                    if compute_Sk
                        if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                            dγk[i1y] = S123 * (n1 + n2) / Det
                            dηk[i1y] = S123 * (n1 * n2) / Det
                            @assert dηk[i1y] >= 0 "2 neg vals dηk[i1y] = $dηk"
                        else
                            F1[i1y] = S123 * (n1 * n2 - n1 * nk[ix, iy] - n2 * nk[ix, iy]) / Det
                        end
                    else
                        if compute_γk
                            dγk[i1y] = S123 * (n1 + n2) / Det
                        end
                        if compute_ηk
                            dηk[i1y] = S123 * (n1 * n2) / Det
                        end
                    end
                end
                if compute_Sk
                    if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                        γk[ix, iy] += integrate_with_log_bins(Run.FSt1[id], iy + 1, i2ky - 1)
                        γk[ix, iy] += integrate_with_log_bins_segment(dγk[i2ky-1], 0.0, kky[i2ky-1], 2 * ky)# add missing part of the integral at the end
                        ηk[ix, iy] += integrate_with_log_bins(Run.FSt2[id], iy + 1, i2ky - 1)
                        ηk[ix, iy] += integrate_with_log_bins_segment(dηk[i2ky-1], 0.0, kky[i2ky-1], 2 * ky)# add missing part of the integral at the end
                        @assert ηk[ix, iy] >= 0 "2 neg vals ηk[ix, iy]"
                    else
                        Sk[ix, iy] += integrate_with_log_bins(Run.FSt[id], iy + 1, i2ky - 1)
                        Sk[ix, iy] += integrate_with_log_bins_segment(F1[i2ky-1], 0.0, kky[i2ky-1], 2 * ky)# add missing part of the integral at the end     
                    end
                else
                    if compute_γk
                        γk[ix, iy] += integrate_with_log_bins(Run.FSt1[id], iy + 1, i2ky - 1)
                        γk[ix, iy] += integrate_with_log_bins_segment(dγk[i2ky-1], 0.0, kky[i2ky-1], 2 * ky)# add missing part of the integral at the end

                    end
                    if compute_ηk
                        ηk[ix, iy] += integrate_with_log_bins(Run.FSt2[id], iy + 1, i2ky - 1)
                        ηk[ix, iy] += integrate_with_log_bins_segment(dηk[i2ky-1], 0.0, kky[i2ky-1], 2 * ky)# add missing part of the integral at the end
                    end
                end

                #Second integral
                F1 .= 0.0
                dγk .= 0.0
                dηk .= 0.0
                #ihky = max(1, 1 - Ik1[1] + ceil(Int, log(ky / 2) / logλ1)) # Mesh point above ky/2
                ihky = ceil(Int, 1+log(ky/(2*kk1[1])) / logλ1) # TODO: Remove previous line if ok
                ihky = clamp(ihky, 1, Mk1)
                # The derminant diverges for k1y->ky/2, but k1x -> \infty. we assume that n(k1x,ky/2)->0 fast enough

                @inbounds for i1y = (ihky):(iy-1) ## n
                    k1y = kky[i1y]
                    k1x = get_k1x_Petviashvilli_Asymp(kx, ky, k1y)
                    k2y = k1y - ky
                    k2x = k1x - kx
                    Det = Det_Petviashvilli_Asymp(ky, k1y)

                    n1 = val_nk(interp_scheeme, Nk, k1x, k1y)
                    n2 = val_nk(interp_scheeme, Nk, k2x, abs(k2y))
                    S123 = S_Petviashvilli(kx, k1x, k2x)
                    if compute_Sk
                        if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                            dγk[i1y] = 2 * S123 * (n2 - n1) / Det
                            dηk[i1y] = 2 * S123 * (n1 * n2) / Det
                            @assert dηk[i1y] >= 0 "3 neg vals dηk[i1y] = $dηk"
                        else
                            F1[i1y] = -2 * S123 * (nk[ix, iy] * n2 - n1 * nk[ix, iy] - n1 * n2) / Det
                        end
                    else
                        if compute_γk
                            dγk[i1y] = 2 * S123 * (n2 - n1) / Det
                        end
                        if compute_ηk
                            dηk[i1y] = 2 * S123 * (n1 * n2) / Det
                        end
                    end
                end
                if compute_Sk
                    if compute_γk || compute_ηk #we check if ned γk or ηk, other we copute St directly
                        γk[ix, iy] += integrate_with_log_bins(Run.FSt1[id], ihky + 1, iy)
                        γk[ix, iy] += integrate_with_log_bins_segment(0.0, dγk[ihky], ky / 2, kky[ihky])# add missing part of the integral at the first point
                        ηk[ix, iy] += integrate_with_log_bins(Run.FSt2[id], ihky + 1, iy)
                        ηk[ix, iy] += integrate_with_log_bins_segment(0.0, dηk[ihky], ky / 2, kky[ihky])# add missing part of the integral at the first point
                        @assert ηk[ix, iy] >= 0 "3 neg vals ηk[ix, iy]"
                    else
                        Sk[ix, iy] += integrate_with_log_bins(Run.FSt[id], ihky + 1, iy)
                        Sk[ix, iy] += integrate_with_log_bins_segment(0.0, F1[ihky], ky / 2, kky[ihky])# add missing part of the integral at the first point
                    end
                else
                    if compute_γk
                        γk[ix, iy] += integrate_with_log_bins(Run.FSt1[id], ihky + 1, iy)
                        γk[ix, iy] += integrate_with_log_bins_segment(0.0, dγk[ihky], ky / 2, kky[ihky])# add missing part of the integral at the first point
                    end
                    if compute_ηk
                        ηk[ix, iy] += integrate_with_log_bins(Run.FSt2[id], ihky + 1, iy)
                        ηk[ix, iy] += integrate_with_log_bins_segment(0.0, dηk[ihky], ky / 2, kky[ihky])# add missing part of the integral at the first point
                    end
                end

                F1[iy] = 0.0 # k2x[iy] = 0
                dγk[iy] = 0.0 # k2x[iy] = 0
                dηk[iy] = 0.0 # k2x[iy] = 0
                @inbounds for i1y = (iy+1):My # In the next integral we make k1y-> -k1y to bring the integral to [ky ,infty]
                    k1y = kky[i1y]
                    k1x = get_k1x_Petviashvilli_Asymp(kx, ky, -k1y)
                    k2y = -k1y - ky
                    k2x = k1x - kx
                    Det = Det_Petviashvilli_Asymp(ky, -k1y)

                    n1 = val_nk(interp_scheeme, Nk, k1x, k1y)
                    n2 = val_nk(interp_scheeme, Nk, k2x, abs(k2y))
                    S123 = S_Petviashvilli(kx, k1x, k2x)
                    if compute_Sk
                        if compute_γk || compute_ηk #we check if ned γk or ηk, other we compute St directly
                            dγk[i1y] = 2 * S123 * (n2 - n1) / Det
                            dηk[i1y] = 2 * S123 * (n1 * n2) / Det
                            @assert dηk[i1y] >= 0 "4 neg vals dηk[i1y] = $dηk"
                        else
                            F1[i1y] = -2 * S123 * (nk[ix, iy] * n2 - n1 * nk[ix, iy] - n1 * n2) / Det
                        end
                    else
                        if compute_γk
                            dγk[i1y] = 2 * S123 * (n2 - n1) / Det
                        end
                        if compute_ηk
                            dηk[i1y] = 2 * S123 * (n1 * n2) / Det
                        end
                    end
                end
                if compute_Sk
                    if compute_γk || compute_ηk #we check if ned γk or ηk, other we compute St directly
                        γk[ix, iy] += integrate_with_log_bins(Run.FSt1[id], iy + 1, My)
                        ηk[ix, iy] += integrate_with_log_bins(Run.FSt2[id], iy + 1, My)
                        @assert ηk[ix, iy] >= 0 "4 neg vals ηk[ix, iy]"
                        γk[ix, iy] *= Sfactor
                        ηk[ix, iy] *= Sfactor
                        Sk[ix, iy] = -γk[ix, iy] * nk[ix, iy] + ηk[ix, iy]
                    else
                        Sk[ix, iy] += integrate_with_log_bins(Run.FSt[id], iy + 1, My)
                        Sk[ix, iy] *= Sfactor
                    end
                else
                    if compute_γk
                        γk[ix, iy] += integrate_with_log_bins(Run.FSt1[id], iy + 1, My)
                        γk[ix, iy] *= Sfactor
                    end
                    if compute_ηk
                        ηk[ix, iy] += integrate_with_log_bins(Run.FSt2[id], iy + 1, My)
                        ηk[ix, iy] *= Sfactor
                    end
                end
            end
        end
    end
    return nothing
end