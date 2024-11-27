# Here we define the matrices and dispersion relation and resonant manifold of Stratified runs


# -------- 3D Stratified_Asymp model --------------
@doc raw"""
    ω_Stratified_Asymp(kh, kz; N=1.0)

Internal gravity wave frequency in the hydrostatic limit
* `kh`: horizontal wave vector modulus
* `kz`: vertical wave vector
* `N=1.0`: buoyancy frequency
"""
function ω_Stratified_Asymp(kh, kz; N=1.0)
    return N * kh / abs(kz)
end

function k1hk2h_vs(kh,p,q)
    k1h = (kh+p+q)/2
    k2h = (kh-p+q)/2
    return k1h, k2h
end

# Lvov-Tabak Formalism
@doc raw"""
    Vk_Strat_Asymp(kh,kz,k1h,k1z,k2h,k2z)

Interaction coefficient on the resonant manifold ``{\bf k} = {\bf k}_1 + {\bf k}_2``, ``\omega_{\bf k} = \omega_1 + \omega_2`` 
* `kh`, `k1h`, `k2h`: horizontal wave vector modulus
* `kz`, `k1z`, `k2z`: vertical wave vector
"""
function Vk_Strat_Asymp(kh,kz,k1h,k1z,k2h,k2z)
    k1hck2h = (kh^2 - k1h^2 - k2h^2)/2
    khck1h = (kh^2 + k1h^2 - k2h^2)/2
    khck2h = (kh^2 - k1h^2 + k2h^2)/2
  
    V = k1hck2h * sqrt( abs( kz/(k1z*k2z) ) ) / (k1h * k2h)
    V = V + khck2h * sqrt( abs( k1z/(kz*k2z) ) ) / (kh * k2h)
    V = V + khck1h * sqrt( abs( k2z/(kz*k1z) ) ) / (kh * k1h)
    V = V * sqrt(kh * k1h * k2h / 32)
    return V
end   

@doc raw"""
    V1_Strat_Asymp(kh,kz,k1h,k1z,k2h,k2z)

Interaction coefficient on the resonant manifold ``{\bf k}_1 = {\bf k} + {\bf k}_2``, ``\omega_1 = \omega_{\bf k} + \omega_2`` 
* `kh`, `k1h`, `k2h`: horizontal wave vector modulus
* `kz`, `k1z`, `k2z`: vertical wave vector
"""
function V1_Strat_Asymp(kh,kz,k1h,k1z,k2h,k2z)
    # Interaction coefficient on the resonant manifold k1 = k + k2 
    k1hck2h = (-kh^2 + k1h^2 + k2h^2)/2
    khck1h = (kh^2 + k1h^2 - k2h^2)/2
    khck2h = (-kh^2 + k1h^2 - k2h^2)/2 

    V = k1hck2h * sqrt( abs( kz/(k1z*k2z) ) ) / (k1h * k2h)
    V = V + khck2h * sqrt( abs( k1z/(kz*k2z) ) ) / (kh * k2h)
    V = V + khck1h * sqrt( abs( k2z/(kz*k1z) ) ) / (kh * k1h)
    V = V * sqrt(kh * k1h * k2h / 32)
    return V
end 

@doc raw"""
    k1zkp_vs(kh,kz,k1h,k2h)

``k_{1z}`` solution of ``{\bf k} = {\bf k}_1 + {\bf k}_2``, ``\omega_{\bf k} = \omega_1 + \omega_2`` (branch +):

``k_{1z} = \frac{k_z}{2k_h} \left( k_h + k_{1h} + k_{2h} + \sqrt{ (k_h + k_{1h} + k_{2h})^2 - 4 k_h k_{1h} } \right).``
* `kh`, `k1h`, `k2h`: horizontal wave vector modulus
* `kz`: vertical wave vector
"""
function k1zkp_vs(kh,kz,k1h,k2h)
    return kz * (kh + k1h + k2h + sqrt( (kh + k1h + k2h)^2 - 4 * kh * k1h)) / (2 * kh) # (* equation 21a of Lvov et al. 2010 *)
end

@doc raw"""
    k1zkm_vs(kh,kz,k1h,k2h)

``k_{1z}`` solution of ``{\bf k} = {\bf k}_1 + {\bf k}_2``, ``\omega_{\bf k} = \omega_1 + \omega_2`` (branch -):

``k_{1z} = \frac{k_z}{2k_h} \left( k_h - k_{1h} - k_{2h} - \sqrt{ (k_h - k_{1h} - k_{2h})^2 + 4 k_h k_{1h} } \right).``
* `kh`, `k1h`, `k2h`: horizontal wave vector modulus
* `kz`: vertical wave vector
"""
function k1zkm_vs(kh,kz,k1h,k2h)
    return kz * (kh - k1h - k2h - sqrt((kh - k1h - k2h)^2 + 4 * kh * k1h)) / (2 * kh) # (* equation 21b of Lvov et al. 2010 *)
end

@doc raw"""
    k1z1p_vs(kh,kz,k1h,k2h)

``k_{1z}`` solution of ``{\bf k}_1 = {\bf k} + {\bf k}_2``, ``\omega_1 = \omega_{\bf k} + \omega_2`` (branch +):

``k_{1z} = \frac{k_z}{2k_h} \left( k_h + k_{1h} + k_{2h} - \sqrt{ (k_h + k_{1h} + k_{2h})^2 - 4 k_h k_{1h} } \right).``
* `kh`, `k1h`, `k2h`: horizontal wave vector modulus
* `kz`: vertical wave vector
"""
function k1z1p_vs(kh,kz,k1h,k2h)
    return kz * (kh + k1h + k2h - sqrt((kh + k1h + k2h)^2 - 4 * kh * k1h)) / (2 * kh) # (* equation 22a of Lvov et al. 2010 *)
end

@doc raw"""
    k1z1m_vs(kh,kz,k1h,k2h)

``k_{1z}`` solution of ``{\bf k}_1 = {\bf k} + {\bf k}_2``, ``\omega_1 = \omega_{\bf k} + \omega_2`` (branch -):

``k_{1z} = \frac{k_z}{2k_h} \left( k_h - k_{1h} + k_{2h} - \sqrt{ (-k_h + k_{1h} - k_{2h})^2 + 4 k_h k_{1h} } \right).``
* `kh`, `k1h`, `k2h`: horizontal wave vector modulus
* `kz`: vertical wave vector
"""
function k1z1m_vs(kh,kz,k1h,k2h)
    return kz * (kh - k1h + k2h - sqrt((-kh + k1h - k2h)^2 + 4 * kh * k1h)) / (2 * kh) # (* equation 22b of Lvov et al. 2010 *)
end


@doc raw"""
    dg_vs(k1h,k1z,k2h,k2z)

Jacobian arizing from resonance in frequency:

``g' = k_{1h} {\rm sign}(k_{1z})/k_{1z}^2 - k_{2h} {\rm sign}(k_{2z})/k_{2z}^2.``
* `k1h`, `k2h`: horizontal wave vector modulus
* `k1z`, `k2z`: vertical wave vectors
"""
function dg_vs(k1h,k1z,k2h,k2z)
    return k1h * sign(k1z) / k1z^2 - k2h * sign(k2z) / k2z^2 # (* Dematteis & Lvov 2021 *)
end



@doc raw"""
    L_vs(ih,iz,p,q,Nk,val_nk,interp_scheeme) 

Integrand of the collision integral multiplied by ``Δ``.
* `ih`, `iz`: external wave vector indices
* `p`, `q`: coordinates on the kinematic box
* `Nk`: wave action spectrum
* `val_nk`, `interp_scheeme`: used for interpolation
"""
function L_vs(ih,iz,p,q,Nk,val_nk,interp_scheeme) 
    # Integrand of the collision integral
    nk = Nk.nk
    kh = Nk.kkh[ih]
    kz = Nk.kkz[iz]
    k1h, k2h = k1hk2h_vs(kh,p,q)
                        
    # k = 1 + 2, branch +
    k1z = k1zkp_vs(kh,kz,k1h,k2h)
    k2z = kz - k1z
    n1 = val_nk(interp_scheeme,Nk,k1h,abs(k1z))
    n2 = val_nk(interp_scheeme,Nk,k2h,abs(k2z)) 
    Vk = Vk_Strat_Asymp(kh,kz,k1h,k1z,k2h,k2z)
    dg = dg_vs(k1h,k1z,k2h,k2z)
    L = Vk^2 * (n1*n2-nk[ih,iz]*(n1+n2))/abs(dg)

    # k = 1 + 2, branch -
    k1z = k1zkm_vs(kh,kz,k1h,k2h)
    k2z = kz - k1z
    n1 = val_nk(interp_scheeme,Nk,k1h,abs(k1z))
    n2 = val_nk(interp_scheeme,Nk,k2h,abs(k2z))
    Vk = Vk_Strat_Asymp(kh,kz,k1h,k1z,k2h,k2z)
    dg = dg_vs(k1h,k1z,k2h,k2z)
    L += Vk^2 * (n1*n2-nk[ih,iz]*(n1+n2))/abs(dg)

    # 1 = k + 2, branch +
    k1z = k1z1p_vs(kh,kz,k1h,k2h)
    k2z = k1z - kz
    n1 = val_nk(interp_scheeme,Nk,k1h,abs(k1z))
    n2 = val_nk(interp_scheeme,Nk,k2h,abs(k2z))
    V1 = V1_Strat_Asymp(kh,kz,k1h,k1z,k2h,k2z) # Warning: keep this order for the arguments (not symmetric by permutation because of the conservation of horizontal momentum)
    dg = dg_vs(k1h,k1z,k2h,k2z)
    L -= 2 * V1^2 * (nk[ih,iz]*n2-n1*(nk[ih,iz]+n2))/abs(dg)    

    # 1 = k + 2, branch -
    k1z = k1z1m_vs(kh,kz,k1h,k2h)
    k2z = k1z - kz
    n1 = val_nk(interp_scheeme,Nk,k1h,abs(k1z))
    n2 = val_nk(interp_scheeme,Nk,k2h,abs(k2z))
    V1 = V1_Strat_Asymp(kh,kz,k1h,k1z,k2h,k2z) # Warning: keep this order for the arguments (not symmetric by permutation because of the conservation of horizontal momentum)
    dg = dg_vs(k1h,k1z,k2h,k2z)
    L -= 2 * V1^2 * (nk[ih,iz]*n2-n1*(nk[ih,iz]+n2))/abs(dg)    
    
    L *= k1h * k2h 
    
    ##### The following lines are used for tests TODO: Remove it or comment it #####
    #L = 1.0
    #L = q
    #L = p
    #L = p^2
    ################################################################################

    return L
end