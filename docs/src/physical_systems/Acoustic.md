## Acoustic2D solver

### Description

In this case, the collisional terms simpilfies to

```math
St_{k}=\frac{2 \pi}{\sqrt{6}a k}\int \left(\mathcal{N}^k_{12}-\mathcal{N}^1_{2k}-\mathcal{N}^2_{k1}\right)\mathrm{d}k_1\mathrm{d}k_2
```

With $\mathcal{N}^k_{12}=|V_{k12}|^2(n_1n_2-n_1n_k-n_kn_2)\delta(k-k_1-k_2)$ and $V_{k12}=V_0\sqrt{kk_1k_2}$.

This is WKE results from assuming a disperse regularisation in the dispersion relation $\omega_k=ck(1+a^2k^2)$, with $ak\ll1$. See [Griffin *et al.*, Phys. Rev. Lett. 128, 224501 (2022)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.128.224501).

### Solver

The computation of the collision integral is relatively simple. We use the resonance condition $k - k_1 - k_2 = 0$ and the symmetry $k_1 \leftrightarrow k_2$ to reduce the collision integral to

```math
St_{k} = I_{k} - 2 J_{k}
```

where

```math
I_{k} = \frac{2 \pi}{\sqrt{6}a k} \int\limits_{0}^{k} ~ |V_{k12}|^2(n_1n_2-n_1n_k-n_kn_2) ~ \mathrm{d}k_1
```

and

```math
J_{k} = \frac{2 \pi}{\sqrt{6}a k} \int\limits_{k}^{k_{\rm max}} ~ |V_{k12}|^2(n_kn_2-n_kn_1-n_1n_2) ~ \mathrm{d}k_1
```

with the integrands being computed using $k_2 = |k - k_1|$ and $k_{\rm max}$ the maximal wave vector modulus.

The wave frequency is $\omega_{k} = ck$ and we use logarithmic grid [`wave_spectrum`](@ref). The different options are described in the [`Acoustic2D`](@ref) documentation.

### Diagnostics

The standard diagnostics should work for this system.

### Theoretical predictions

The Kolmogorov-Zakharov energy spectrum is 

```math
e_k = 2\pi k \omega_k n_{\bf k} = \frac{4^{3/2} \sqrt{a c}}{3 \pi} ~ k^{-1}.
```


## Acoustic3Dsolver

### Description

In this case, the collisional terms simpilfies to

```math
St_{k}=\frac{4 \pi^2}{c k^2}\int \left(\mathcal{N}^k_{12}-\mathcal{N}^1_{2k}-\mathcal{N}^2_{k1}\right)k_1 k_2 \mathrm{d}k_1\mathrm{d}k_2
```

With $\mathcal{N}^k_{12}=|V_{k12}|^2(n_1n_2-n_1n_k-n_kn_2)\delta(k-k_1-k_2)$ and $V_{k12}=V_0\sqrt{kk_1k_2}$.

This WKE describes 3D acoustic waves with frequency $\omega_k=ck$. See [Zakharov and Sagdeev, Dokl. Akad. Nauk SSSR 192, 297--300 (1970)](http://mi.mathnet.ru/dan35409) and [Kochurin and Kuznetsov, arXiv, 2407.08352, (2024)](https://arxiv.org/abs/2407.08352).

### Solver

The computation of the collision integral is the same than for [Acoustic2D solver](@ref), up to the prefactor. We use logarithmic grid [`wave_spectrum`](@ref). The different options are described in the [`Acoustic3D`](@ref) documentation.

### Diagnostics

The standard diagnostics should work for this system.

### Theoretical predictions

The Kolmogorov-Zakharov energy spectrum is the Zakharov-Sagdeev (ZS) spectrum.

```math
e_k = 4\pi k^2 \omega_k n_{\bf k} = \frac{\sqrt{c P} k^{-3/2}}{\sqrt{3\pi(\pi + 4 \ln 2 − 1)}}
```

with $P$ the energy injection rate.






### List of structures and methods for Acoustic solvers
```@autodocs
Modules = [WavKinS]
Order   = [:type, :function]
Pages   = ["src/physical_systems/Acoustic/basics.jl", "src/physical_systems/Acoustic/structure.jl", "src/physical_systems/Acoustic/collision_integral.jl"]
```
