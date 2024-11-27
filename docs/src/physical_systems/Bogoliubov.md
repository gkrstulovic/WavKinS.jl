## Bogoliubov3Dsolver

### Description

In this case, the collisional terms simpilfies to

$$St_{k}= \frac{4 \pi}{ck^2} \int \left(\mathcal{N}^k_{12}-\mathcal{N}^1_{2k}-\mathcal{N}^2_{k1}\right)k_1 k_2 \mathrm{d}k_1\mathrm{d}k_2$$

With $\mathcal{N}^k_{12}=|V_{12}^k|^2(n_1n_2-n_1n_k-n_kn_2)\delta(\omega_k-\omega_1-\omega_2)$, 

$$V_{12}^k = V_0\sqrt{kk_1k_2} \left[ \frac{1}{2\sqrt{\eta_k \eta_1 \eta_2}} + \frac{\sqrt{\eta_k \eta_1 \eta_2}}{6 k k_1 k_2} \left( \frac{k^3}{\eta_k} - \frac{k^3}{\eta_1} - \frac{k^3}{\eta_2} \right) \right],$$

and $\eta_k = \sqrt{1+(k\xi)^2/2}$. $c$ is the speed of sound and $\xi$ the healing length.

The dispersion relation is $\omega_k=ck\sqrt{1+(k\xi)^2/2}$. See [Zhu *et al.*, arXiv, 2408.15163 (2024)](https://arxiv.org/abs/2408.15163). Note that the limit $\xi \rightarrow 0$ corresponds to [Acoustic3Dsolver](@ref).

### Solver

We use the resonance condition $\omega_k - \omega_1 - \omega_2 = 0$ and the symmetry $k_1 \leftrightarrow k_2$ to reduce the collision integral to

$$St_{k} = I_{k} - 2 J_{k}$$

where

$$I_{k} = \frac{4 \pi}{ck^2} \int\limits_{0}^{k} ~ \frac{|V_{k12}|^2}{g(k_2)} (n_1n_2-n_1n_k-n_kn_2) ~ k_1 k_2 ~ \mathrm{d}k_1$$

and

$$J_{k} = \frac{4 \pi}{ck^2} \int\limits_{k}^{k_{\rm max}} ~ \frac{|V_{k12}|^2}{g(k_2)}(n_kn_2-n_kn_1-n_1n_2) ~ k_1 k_2 ~ \mathrm{d}k_1$$

with the integrands being computed using $k_2$ solution to $\omega_2 = \pm (\omega_k - \omega_1)$, i.e.

$$k_2 = k_2(k, k_1) = \frac{1}{\xi} \left[ \sqrt{1 + \frac{2 \xi^2}{c^2} (\omega_k - \omega_1)^2} - 1 \right]^{1/2},$$

and $k_{\rm max}$ the maximal wave vector modulus. We introduced the factor 

$$g(k_2) = \left| \frac{\partial (\omega_k - \omega_1 - \omega_2)}{\partial k_2} \right| = \left| \frac{\partial (\omega_1 - \omega_k - \omega_2)}{\partial k_2} \right| = c
\frac{1 + \xi^2 k^2}{\sqrt{1 + \xi^2 k^2/2}}$$

comming from the integration of  $\delta(\omega_k - \omega_2 - \omega_2)$.

We use logarithmic grid [`wave_spectrum`](@ref). The different options are described in the [`Bogoliubov3D`](@ref) documentation.

### Diagnostics

The standard diagnostics should work for this system.

### Theoretical predictions

For large scale ($k\xi \ll 1$), the system describes acoustic waves and the Kolmogorov-Zakharov energy spectrum is the Zakharov-Sagdeev (ZS) spectrum.

$$E(k) = 4\pi k^2 \omega_k n_{\bf k} = \frac{\sqrt{c P} k^{-3/2}}{\sqrt{3\pi(\pi + 4 \ln 2 − 1)}}$$

with $P$ the energy injection rate.

For small scale ($k\xi \gg 1$), the Kolmogorov-Zakharov energy spectrum is 

$$E(k) = 4\pi k^2 \omega_k n_{\bf k} = \frac{ 2^{3/4} \sqrt{c P \xi^5} k}{\sqrt{\pi(\pi − 4 \ln 2)}}.$$



 
 
 
 
 
 
 
 
  


 


### List of structures and methods for Bogoliubov solvers
```@autodocs
Modules = [WavKinS]
Order   = [:type, :function]
Pages   = ["src/physical_systems/Bogoliubov/basics.jl", "src/physical_systems/Bogoliubov/structure.jl", "src/physical_systems/Bogoliubov/collision_integral.jl"]
```