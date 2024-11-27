using NCDatasets
#using DataStructures


@doc raw"""
    init_IO(Run, outputDir::String)

Inputs/Outputs initialization.
* `Run`: WavKinS simulation structure
* `outputDir`: path of the output directory 
"""
function init_IO(Run, outputDir::String)
    if Run.Nk_arguments == 1
        IsotropicRun = true
    else
        IsotropicRun = false
    end
    NCDataset(outputDir * "WKE_" * Run.name * "_data.nc", "c") do ds

        ds.attrib["Dataset"] = Run.name
        ds.attrib["Created with git commit id"] = commit_version

        dsGlobal = defGroup(ds, "Global")
        dsSP = defGroup(ds, "Spectral")
        if IsotropicRun
            defDim(ds, "k", length(Run.Nk.kk))
        else
            defDim(ds, "kh", Run.Nk.Mh)
            defDim(ds, "kz", Run.Nk.Mz)
        end
        defDim(ds, "timesSP", Inf)
        defDim(ds, "times", Inf)

        if IsotropicRun
            k_var = defVar(dsSP, "k", Float64, ("k",))
            k_var.attrib["long_name"] = "wave vectors"
            k_var[:] = Run.Nk.kk
        else
            kh_var = defVar(dsSP, "kh", Float64, ("kh",))
            kh_var.attrib["long_name"] = "kh wave vectors"
            kh_var[:] = Run.Nk.kkh

            kz_var = defVar(dsSP, "kz", Float64, ("kz",))
            kz_var.attrib["long_name"] = "kz wave vectors"
            kz_var[:] = Run.Nk.kkz
        end

        # N_var = defVar(dsGlobal, "N", Float64, ("times",))
        # N_var.attrib["long_name"] = "wave action"
        for (varname, glob_output) in Run.diags.glob_diag
            glob_var = defVar(dsGlobal, varname, Float64, ("times",))
            glob_var.attrib["long_name"] = glob_output.longname
        end


        timesSP_var = defVar(dsSP, "TimesSP", Float64, ("timesSP",))
        timesSP_var.attrib["long_name"] = "Times of spectral outputs"

        if IsotropicRun
            # nk_var = defVar(dsSP, "nk", Float64, ("k", "timesSP"))
            # nk_var.attrib["long_name"] = "wave action spectrum"
            for (varname, sp_outs) in Run.diags.sp_outs
                sp_var = defVar(dsSP, varname, Float64, ("k", "timesSP"))
                sp_var.attrib["long_name"] = sp_outs.longname
            end

        else
            # nk_var = defVar(dsSP, "nk", Float64, ("kh", "kz", "timesSP"))
            # nk_var.attrib["long_name"] = "wave action spectrum"
            for (varname, sp_outs) in Run.diags.sp_outs
                sp_var = defVar(dsSP, varname, Float64, ("kh", "kz", "timesSP"))
                sp_var.attrib["long_name"] = sp_outs.longname
            end
        end

    end
end


@doc raw"""
    output_global(Run, outputDir::String)

Write global output (total wave action, total energy, ...).
* `Run`: WavKinS simulation structure
* `outputDir`: path of the output directory  
"""
function output_global(Run, outputDir::String)
    NCDataset(outputDir * "WKE_" * Run.name * "_data.nc", "a") do ds
        nt = ds.dim["times"]
        dsGlobal = ds.group["Global"]
        for (varname, glob_output) in Run.diags.glob_diag
            dsGlobal[varname][nt+1] = glob_output.out[end]
        end

    end

end

@doc raw"""
    output_spectra(Run, outputDir::String)

Write spectra.
* `Run`: WavKinS simulation structure
* `outputDir`: path of the output directory   
"""
function output_spectra(Run, outputDir::String)
    NCDataset(outputDir * "WKE_" * Run.name * "_data.nc", "a") do ds
        nt = ds.dim["timesSP"]
        dsSP = ds.group["Spectral"]

        for (varname, sp_outs) in Run.diags.sp_outs
            if sp_outs.write_sp
                if Run.Nk_arguments == 1
                    dsSP[varname][:, nt+1] = sp_outs.sp
                elseif Run.Nk_arguments == 2
                    dsSP[varname][:, :, nt+1] = sp_outs.sp
                else
                    @error "Physical systems with more than 2 degree of freedom are not implemented"
                end
                dsSP["TimesSP"][nt+1] = Run.t
            end
        end
    end
end

@doc raw"""
    load_spectrum(Run, outputDir::String; irestart=-1)

Load wave action spectrum.
* `Run`: WavKinS simulation structure
* `outputDir`: path of the output directory   
* `irestart`: time step
By default, it loads the spectrum at the latest time step.
"""
function load_spectrum(Run, outputDir::String; irestart=-1)

    NCDataset(outputDir * "WKE_" * Run.name * "_data.nc", "r") do ds
        dsSP = ds.group["Spectral"]
        nk_var = dsSP["nk"]
        if irestart < 0
            irestart = ds.dim["timesSP"]
        end
        timesSP_var = dsSP["TimesSP"]
        if Run.Nk_arguments == 1
            @. Run.Nk.nk = nk_var[:, irestart]
        elseif Run.Nk_arguments == 2 
            @. Run.Nk.nk = nk_var[:,:, irestart]
        else
            return error
        end
        Run.t = timesSP_var[irestart]
    end
    return nothing
end

@doc raw"""
    load_spectral_for_change_mesh!(Nk::wave_spectrum, datafile::String, Run; irestart=-1)

Load wave action spectrum to change the mesh.
* `Run`: WavKinS simulation structure
* `datafile`: path of the output file  
* `irestart`: time step
By default, it loads the spectrum at the latest time step.
"""
function load_spectral_for_change_mesh!(Nk::wave_spectrum, datafile::String, Run; irestart=-1)

    NCDataset(datafile, "r") do ds
        dsSP = ds.group["Spectral"]
        nk_var = dsSP["nk"]
        if irestart < 0
            irestart = ds.dim["timesSP"]
        end
        timesSP_var = dsSP["TimesSP"]
        @. Nk.nk = nk_var[:, irestart]
        Run.t = timesSP_var[irestart]
    end
    return nothing
end
