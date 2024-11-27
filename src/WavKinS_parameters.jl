@doc raw"""
Basic structure containing information about the simulation.

    dt::Float64 # time step of the simulation
    tplot::Float64 # Plot every tplot times
    tglobal::Float64 # Compute, store and write global quantities every tglobal times 
    tspstore::Float64 # Compute and store spectral quantities every tspstore times
    tspwrite::Float64 # write spectra every tspwrite

    outputDir::String #output directory
    write_global::Bool #write global
    write_spectral::Bool #write spectra
    
"""
mutable struct simulation_parameters
    dt::Float64 # time step of the simulation
    tplot::Float64 # Plot every tplot times
    tglobal::Float64 # Compute, store and write global quantities every tglobal times 
    tspstore::Float64 # Compute and store spectral quantities every tspstore times
    tspwrite::Float64 # write spectra every tspwrite
    
    outputDir::String #output directory
    write_global::Bool #write global
    write_spectral::Bool #write spectra
end

@doc raw"""
    simulation_parameters(dt,tplot,tglobal,tspstore)

Constructor for [`simulation_parameters`](@ref "simulation_parameters") structure.
* `dt`: time step of the simulation
* `tplot`: plot every tplot times
* `tglobal`: compute, store and write global quantities every tglobal times 
* `tspstore`: compute and store spectral quantities every tspstore times
"""
function simulation_parameters(dt,tplot,tglobal,tspstore)
    simulation_parameters(dt,tplot,tglobal,tspstore,Inf,"./",false,false)
end




#TODO: I think tspstore can be removed