using Base.Threads
# Here we define collisional integrals for Stratified runs

@doc raw"""
    St_k!(Run::Stratified_Asymp)

Compute collision integral for [`Stratified_Asymp`](@ref "Stratified_Asymp").

"""
function St_k!(Run::Stratified_Asymp)
    Nk = Run.Nk
    kkh = Nk.kkh
    kkz = Nk.kkz
    Sk = Run.Sk.nk

    interp_scheeme = Run.interp_scheeme
    update_coeff_interp!(interp_scheeme, Nk)

    Sfactor = 4 * pi # Take into account the jacobian for (k1h,k2h) → (p,q)

    @sync for id = 1:Threads.nthreads()
        inds = Run.partition[id]

        Fp = Run.FStp[id]
        Fq = Run.FStq[id]
        Mq = Run.kin_box.Mq
        qq = Run.kin_box.qq
        λq = Run.kin_box.λq
        logλq = Run.kin_box.logλq

        @spawn begin
            @inbounds for iz in inds[2], ih in inds[1]

                Sk[ih,iz] = 0.0
                kh = kkh[ih]
                kz = kkz[iz]
                aa = Run.kin_box.aa[ih]
                Ma = Run.kin_box.Ma[ih]
                λa = Run.kin_box.λa[ih]
                logλa = Run.kin_box.logλa[ih]

                χ = aa[1]/kh
                δ = qq[1]/kh

                # Part p ∈ ]-kh:0]
                @inbounds for ia in eachindex(aa)
                    p = aa[ia]-kh
                    for iq in 1:Mq    
                        q = qq[iq]
                        Fq[iq] = L_vs(ih,iz,p,q,Nk,val_nk,interp_scheeme) / Δ_vs(kkh[ih],p,q)
                    end
                    Fp[ia] = integrate_with_log_bins(Mq,λq,logλq,qq,Fq,2,Mq) # imin = 2 to not compute q∈[0:qmin] cells by interpolation 
                    # Singularities p≠±kh q=0 
                    Fp[ia] += (L_vs(ih,iz,p,0,Nk,val_nk,interp_scheeme) + L_vs(ih,iz,p,qq[1],Nk,val_nk,interp_scheeme)) * sqrt(2*δ/(kh^2-p^2)) 
                end
                Sk[ih,iz] = integrate_with_log_bins(Ma,λa,logλa,aa,Fp,2,Ma) # imin = 2 to not compute p∈[-kh:-kh+amin] cells by interpolation 

                # Part p ∈ [0:kh[
                @inbounds for ia in 1:Ma-1 # We end at Ma-1 because ia=length(aa) is p=0 (computed in the previous part)
                    p = kh-aa[ia]
                    for iq in 1:Mq    
                        q = qq[iq]
                        Fq[iq] = L_vs(ih,iz,p,q,Nk,val_nk,interp_scheeme) / Δ_vs(kkh[ih],p,q)
                    end
                    Fp[ia] = integrate_with_log_bins(Mq,λq,logλq,qq,Fq,2,Mq) # imin = 2 to not compute q∈[0:qmin] cells by interpolation 
                    # Singularities p≠±kh q=0 
                    Fp[ia] += (L_vs(ih,iz,p,0,Nk,val_nk,interp_scheeme) + L_vs(ih,iz,p,qq[1],Nk,val_nk,interp_scheeme)) * sqrt(2*δ/(kh^2-p^2)) 
                end
                Sk[ih,iz] += integrate_with_log_bins(Ma,λa,logλa,aa,Fp,2,Ma) # imin = 2 to not compute p∈[kh-amin:kh] cells by interpolation 

                # Singularities p=±kh q≠0 
                for iq in 1:Mq
                    q = qq[iq]
                    Fq[iq] = L_vs(ih,iz,-kh,q,Nk,val_nk,interp_scheeme) + L_vs(ih,iz,-kh*(1-χ),q,Nk,val_nk,interp_scheeme)
                    Fq[iq] += L_vs(ih,iz,kh*(1-χ),q,Nk,val_nk,interp_scheeme) + L_vs(ih,iz,kh,q,Nk,val_nk,interp_scheeme)
                    Fq[iq] *= sqrt(2*χ/(q*(2*kh+q)))
                end
                Sk[ih,iz] += integrate_with_log_bins(Mq,λq,logλq,qq,Fq,2,Mq)
    
                # Singularities p=±kh q=0 here
                Sk[ih,iz] += (L_vs(ih,iz,kh*(-1+χ),qq[1],Nk,val_nk,interp_scheeme) + L_vs(ih,iz,kh*(1-χ),qq[1],Nk,val_nk,interp_scheeme)) * 2 * (pi - 2*asin(1-χ)) * asin(sqrt(δ/2)) 
                
                # Collisional integral prefactor
                Sk[ih,iz] *= Sfactor  
                                 
            end
        end
    end
    return nothing
end
