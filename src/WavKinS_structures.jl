## Grid and spectra

@doc raw"""
Basic structure for a 1D field defined on an arbitrary grid. 
    
    M::Int: number of grid points
    kk::Vector{Float64}: grid
    F::Vector{Float64}: field

"""
struct field_grid_1D
    M::Int # number of points
    kk::Vector{Float64} # grid
    F::Vector{Float64} # field
end

@doc raw"""
    field_grid_1D(kk::Vector{Float64})

Constructor of the [`field_grid_1D`](@ref "field_grid_1D") structure. 
* `kk`: grid
"""
function field_grid_1D(kk::Vector{Float64})
    M = length(kk)
    F = zeros(M)
    field_grid_1D(M,kk,F)
end


@doc raw"""
Basic structure for an isotropic or 1D field with logarithmic grid. 
    
    M::Int # number of points
    nk::Vector{Float64} # wave action
    λ::Float64 # logarithmic increment of the mesh
    logλ::Float64 # log(λ)
    kk::Vector{Float64} # wave vector

"""
struct wave_spectrum
    M::Int # number of points
    nk::Vector{Float64} # wave action
    λ::Float64 # logarithmic increment of the mesh
    logλ::Float64 # log(λ)
    kk::Vector{Float64} # wave vector
end

@doc raw"""
    wave_spectrum(kmin::Float64, kmax::Float64, M::Int)

Constructor of a [`wave_spectrum`](@ref "wave_spectrum") structure. 
* `kmin`: minimal wave vector 
* `kmax`: maximal wave vector
* `M`: number of grid points
The wave vector grid points are ``k[i] = k_{\rm min} λ^{i-1}`` where ``λ`` is the logarithmic increment such that ``k[M] = k_{\rm max}``.
"""
function wave_spectrum(kmin::Float64, kmax::Float64, M::Int)
    kk = LogRange(kmin,kmax,M)
    λ = kk[2]/kk[1]
    nk = zeros(M)
    wave_spectrum(M, nk, λ, log(λ), kk)
end

@doc raw"""
Basic structure for an axisymmetric or 2D (with ``k_z ↔ -k_z`` symmetry) spectrum. 
    
    Mh::Int # number of points in the horizontal
    Mz::Int # number of points in the vertical 
    nk::Array{Float64,2} # wave action
    λh::Float64 # logarithmic increment of the mesh in the horizontal
    λz::Float64 # logarithmic increment of the mesh in the vertical
    logλh::Float64 # log(λh)
    logλz::Float64 # log(λz)
    kkh::Vector{Float64} # horizontal wave vector modulus
    kkz::Vector{Float64} # vertical wave vector modulus
    kk::Array{Float64,2} # wave vector modulus

"""
struct wave_spectrum_khkz
    Mh::Int # number of points in the horizontal
    Mz::Int # number of points in the vertical 
    nk::Array{Float64,2} # wave action
    λh::Float64 # logarithmic increment of the mesh in the horizontal
    λz::Float64 # logarithmic increment of the mesh in the vertical
    logλh::Float64 # log(λh)
    logλz::Float64 # log(λz)
    kkh::Vector{Float64} # horizontal wave vector modulus
    kkz::Vector{Float64} # vertical wave vector modulus
    kk::Array{Float64,2} # wave vector modulus
end

@doc raw"""
    wave_spectrum_khkz(khmin::Float64, khmax::Float64, Mh::Int, kzmin::Float64, kzmax::Float64, Mz::Int)

Constructor of a [`wave_spectrum_khkz`](@ref "wave_spectrum_khkz") structure. 
* `khmin`: minimal horizontal wave vector 
* `khmax`: maximal horizontal wave vector
* `Mh`: number of horizontal wave vector grid points
* `kzmin`: minimal vertical wave vector 
* `kzmax`: maximal vertical wave vector
* `Mz`: number of vertical wave vector grid points
The horizontal wave vector grid points are ``k_h[i_h] = k_{h {\rm min}} λ_h^{i_h-1}`` where ``λ_h`` is the logarithmic increment such that ``k_h[M_h] = k_{h {\rm max}}``.
The vertical wave vector grid points are ``k_z[i_z] = k_{z {\rm min}} λ_z^{i_z-1}`` where ``λ_z`` is the logarithmic increment such that ``k_z[M_z] = k_{z {\rm max}}``.
"""
function wave_spectrum_khkz(khmin::Float64, khmax::Float64, Mh::Int, kzmin::Float64, kzmax::Float64, Mz::Int)
    kkh = LogRange(khmin,khmax,Mh)
    λh = kkh[2]/kkh[1]
    kkz = LogRange(kzmin,kzmax,Mz)
    λz = kkz[2]/kkz[1]
    nk = zeros(Float64, Mh, Mz)
    KKH = kkh .* ones(length(kkz))'
    KKZ = ones(length(kkh)) .* kkz'
    kk = (KKH.^2 + KKZ.^2).^0.5
    wave_spectrum_khkz(Mh, Mz, nk, λh, λz, log(λh), log(λz), kkh, kkz, kk)
end

@doc raw"""
Basic structure for a kinematic box with ``p`` and ``q`` variables. 
![](assets/kinematic_box.png)
We use ``a`` such that ``p=a-k_h ∈ [-k_h:0]`` or ``p=k_h-a ∈ [0:k_h]``. 
Note: The horizontal wave vectors of first and second waves are ``k_{1h} = (k_h + p + q)/2`` and ``k_{2h} = (k_h - p + q)/2``.

    Ma::Vector{Int} # number of points in a (depends on kh)
    Mq::Int # number of points in q 
    λa::Vector{Float64} # logarithmic increment of the meshes in a (depends on kh)
    λq::Float64 # logarithmic increment of the mesh in q
    logλa::Vector{Float64} # log(λa)
    logλq::Float64 # log(λq)
    aa::Vector{Vector{Float64}} # meshes of a (depends on kh)
    qq::Vector{Float64} # mesh of q

"""
struct kinematic_box
    Ma::Vector{Int} # number of points in a (depends on kh)
    Mq::Int # number of points in q 
    λa::Vector{Float64} # logarithmic increment of the meshes in a
    λq::Float64 # logarithmic increment of the mesh in q
    logλa::Vector{Float64} # log(λa)
    logλq::Float64 # log(λq)
    aa::Vector{Vector{Float64}} # meshes of a
    qq::Vector{Float64} # mesh of q
end

@doc raw"""
    kinematic_box(amin::Vector{Float64}, amax::Vector{Float64}, Ma::Vector{Int}, qmin::Float64, qmax::Float64, Mq::Int)

Constructor of a [`kinematic_box`](@ref "kinematic_box") structure. 
* `amin`: minimal value of ``a`` (depends on ``k_h``)
* `amax`: maximal value of ``a`` (depends on ``k_h``)
* `Ma`: number of grid points in a (depends on ``k_h``)
* `qmin`: minimal value of ``q``
* `qmax`: maximal value of ``q``
* `Mq`: number of grid points in ``q``
The ``a`` grid points are ``a[i] = a_{\rm min} λ_a^{i-1}`` where ``λ_a`` is the logarithmic increment such that ``a[M_a] = a_{\rm max}``.
The ``q`` grid points are ``q[i] = q_{\rm min} λ_q^{i-1}`` where ``λ_q`` is the logarithmic increment such that ``q[M_q] = q_{\rm max}``.
"""
function kinematic_box(amin::Vector{Float64}, amax::Vector{Float64}, Ma::Vector{Int}, qmin::Float64, qmax::Float64, Mq::Int)
    @assert length(amin) == length(amax) && length(amin) == length(Ma)
    Mh=length(Ma)
    λa=Vector{Float64}(undef,Mh)
    aa=Vector{Vector{Float64}}(undef,Mh)

    for ih in eachindex(Ma)
        aa[ih]=LogRange(amin[ih],amax[ih],Ma[ih])
        λa[ih] = aa[ih][2]/aa[ih][1]
    end
    
    qq=LogRange(qmin,qmax,Mq) 
    λq=qq[2]/qq[1]

    kinematic_box(Ma, Mq, λa, λq, log.(λa), log(λq), aa, qq)
end


## Forcing and dissipation coefficients

@doc raw"""
Structure for storing forcing and dissipation coefficients ``f_{\bf k}`` and ``d_{\bf k}``.

    f::VecOrMat{Float64} # forcing of the WKE
    D::VecOrMat{Float64} # dissipation of the WKE 

"""
mutable struct force_dissipation
    f::VecOrMat{Float64} # forcing of the WKE
    D::VecOrMat{Float64} # dissipation of the WKE 
end

@doc raw"""
    force_dissipation(M::Int)

Constructor of [`force_dissipation`](@ref "force_dissipation") for isotropic or 1D systems.
* `M`: number of wave vector grid points 
"""
function force_dissipation(M::Int)
    force_dissipation(zeros(M), zeros(M))
end

@doc raw"""
    force_dissipation(Mh::Int, Mz::Int)

Constructor of [`force_dissipation`](@ref "force_dissipation") for axisymmetric or 2D systems.
* `Mh`: number of horizontal wave vector grid points 
* `Mz`: number of vertical wave vector grid points
"""
function force_dissipation(Mh::Int, Mz::Int)
    force_dissipation(zeros(Mh, Mz), zeros(Mh, Mz))
end



## Outputs

@doc raw"""
Basic structure for storing a global variable. 
    
    longname::string  # long name
    out::Vector{Float64} # vector containing the output

"""
mutable struct global_ouput
    longname::String  # long name of the variable
    out::Vector{Float64} # vector containing the output
end

@doc raw"""
    global_ouput(longname)

Constructor of a [`global_output`](@ref "global_ouput"). 
* `longname`: long name of the global variable
"""
function global_ouput(longname)
    global_ouput(longname, Array{Float64}(undef, 0))
end


@doc raw"""
Basic structure for storing a spectrum. 
    
    longname::String  # long name of the variable
    sp::VecOrMat{Float64}  # spectrum
    write_sp::Bool # write and compute this spectrum

"""
mutable struct spectral_output
    longname::String  # long name of the variable
    sp::VecOrMat{Float64}  # spectrum
    write_sp::Bool  # write and compute this spectrum
end

@doc raw"""
     spectral_output(longname,M)

Constructor of a [`spectral_output`](@ref "spectral_output") for isotropic or 1D systems.
* `longname`: long name of the spectrum
* `M`: number of wave vector grid points 
"""
function spectral_output(longname::String, M::Int)
    spectral_output(longname, zeros(M), true)
end

@doc raw"""
     spectral_output(longname,Mh, Mz)

Constructor of a [`spectral_output`](@ref "spectral_output") for axisymmetric or 2D systems.
* `longname`: long name of the spectrum
* `Mh`: number of horizontal wave vector grid points 
* `Mz`: number of vertical wave vector grid points
"""
function spectral_output(longname::String, Mh::Int, Mz::Int)
    spectral_output(longname, zeros(Mh, Mz), true)
end


## Diagnostics

@doc raw"""
Structure containing diagnostics, spectra, etc.
    
    glob_diag::Dict{String, global_ouput} # dictionary with global outputs
    sp_outs::Dict{String, spectral_output} # dictionary with spectral outputs
    sp_store::Dict{String, Array{Float64}} # dictionary with stored in memory spectral quantities

"""
mutable struct diagnostic_container
    glob_diag::Dict{String,global_ouput} # dictionary with global outputs
    sp_outs::Dict{String,spectral_output} # dictionary with spectral outputs
    sp_store::Dict{String,Array{Float64}} # dictionary with stored in memory spectral quantities
end

@doc raw"""
    function diagnostic_container(M::Int)

Constructor of a [`diagnostic_container`](@ref "diagnostic_container") for a [`wave_spectrum`](@ref "wave_spectrum") with `M` points.

    The dictionaries are initialised by default with standard outputs:
    glob_diag : It contains the keys "Times", "H", "N", and "Disp":
        "Times" => "Times of global quantities"
        "H" => "Total energy"  
        "N" => "Total wave action"  
        "Disp" => ""Total dissipation""    

    sp_outs : It contains the keys "nk", "Ek", and "Pk"
        "nk" => "Wave action spectrum"
        "Ek" => "Energy spectrum"
        "Pk" => "Energy flux"

    sp_store : It contains the keys "nk", "Pk", and "Tsp"
        "nk" => "Wave action spectrum"
        "Pk" => "Energy spectrum"
        "Tsp" => "Times of stored spectra"
    
"""
function diagnostic_container(M::Int)
    glob_diag = Dict{String,global_ouput}()
    glob_diag["Times"] = global_ouput("Times of global quantities")
    glob_diag["H"] = global_ouput("Total energy")
    glob_diag["N"] = global_ouput("Total wave action")
    glob_diag["Disp"] = global_ouput("Total energy dissipation")

    sp_outs = Dict{String,spectral_output}()
    sp_outs["nk"] = spectral_output("Wave action spectrum", M)
    sp_outs["Ek"] = spectral_output("Energy spectrum", M)
    sp_outs["Pk"] = spectral_output("Energy flux", M)

    sp_store = Dict{String,Array{Float64}}()
    sp_store["nk"] = Array{Float64}(undef, M, 0)
    sp_store["Pk"] = Array{Float64}(undef, M, 0)
    sp_store["Tsp"] = Array{Float64}(undef, 0)

    diagnostic_container(glob_diag, sp_outs, sp_store)
end

@doc raw"""
    function diagnostic_container(Mh::Int, Mz::Int)

Constructor of a [`diagnostic_container`](@ref "diagnostic_container") for [`wave_spectrum_khkz`](@ref "wave_spectrum_khkz") with `Mh` and `Mz` points. 
Same ininitialisation as [`diagnostic_container(M::Int)`](@ref "diagnostic_container(Mh::Int)").

"""
function diagnostic_container(Mh::Int, Mz::Int)
    glob_diag = Dict{String,global_ouput}()
    glob_diag["Times"] = global_ouput("Times of global quantities")
    glob_diag["H"] = global_ouput("Total energy")
    glob_diag["N"] = global_ouput("Total wave action")
    glob_diag["Disp"] = global_ouput("Total dissipation")

    sp_outs = Dict{String,spectral_output}()
    sp_outs["nk"] = spectral_output("Wave action spectrum", Mh, Mz)
    sp_outs["Ek"] = spectral_output("Energy spectrum", Mh, Mz)
    sp_outs["Pkh"] = spectral_output("Energy flux along kh", Mh, Mz)
    sp_outs["Pkz"] = spectral_output("Energy flux along kz", Mh, Mz)

    sp_store = Dict{String,Array{Float64}}()
    sp_store["nk"] = Array{Float64}(undef, Mh, Mz, 0)
    sp_store["Pkh"] = Array{Float64}(undef, Mh, Mz, 0)
    sp_store["Pkz"] = Array{Float64}(undef, Mh, Mz, 0)
    sp_store["Tsp"] = Array{Float64}(undef, 0)

    diagnostic_container(glob_diag, sp_outs, sp_store)
end