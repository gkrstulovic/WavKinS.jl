@doc raw"""
    meshgrid(x, y)

Create 2D grids from the 1D vectors `x` and `y`.

"""
function meshgrid(x, y)
    X = [i for i in x, j in 1:length(y)]
    Y = [j for i in 1:length(x), j in y]
    return X, Y
end

@doc raw"""
    LogRange(kmin,kmax,M)

Create `M` logarithmic grid points between `kmin` and `kmax`.

"""
function LogRange(kmin,kmax,M)
    logkmin = log(kmin)
    logkmax = log(kmax)
    logk = LinRange(logkmin,logkmax,M)
    k = exp.(logk)
    return k
end

@doc raw"""
    change_mesh!(Nknew,Nkold;interp_scheeme=BS_interp)

Fill `Nknew` wave action spectrum from interpolation `interp_scheeme` of `Nkold`.

"""
function change_mesh!(Nknew,Nkold;interp_scheeme=BS_interp)
    interp_scheeme = interp_scheeme(Nkold.kk)
    update_coeff_interp!(interp_scheeme,Nkold)
    for i in eachindex(Nknew.nk)
        Nknew.nk[i] = val_nk(interp_scheeme,Nkold, Nknew.kk[i])
    end
    return nothing
end

@doc raw"""
    area_ratio(yi,yu,yl,yr)

Area of a box (x,y) ∈ [xl,xr] × [`yi`,`yu`] that is below the line passing by (xl,`yl`) and (xr,`yr`), divided by the total area of the box.

Note: does not depend on xl and xr.

"""
function area_ratio(yi,yu,yl,yr)
    if yl <= yi
        if yr <= yi
            R=0.0;
        elseif yr <= yu
            R=(yr-yi).^2 ./ (2*(yr-yl).*(yu-yi));
        else
            R=(yr-(yi+yu)/2)./(yr-yl);
        end
    elseif yl <= yu
        if yr <= yi
            R=(yl-yi).^2 ./ (2*(yl-yr).*(yu-yi));
        elseif yr <= yu
            R=((yl+yr)/2-yi) ./ (yu-yi);
        else
            R=1-(yu-yl).^2 ./ (2*(yr-yl).*(yu-yi));
        end
    else
        if yr <= yi
            R=(yi+yu-2*yl) ./ (2*(yr-yl));
        elseif yr <= yu
            R=1-(yu-yr).^2 ./ (2*(yl-yr).*(yu-yi));
        else
            R=1.0;
        end
    end
    return R
end

@doc raw"""
    area_ratio_grid(x, y, f)

For each cell of the grid `x` × `y`, return the area under `f` divided by the total area of the box. 

Note: convolution allows to integrate over `y` < `f`.

"""
function area_ratio_grid(x, y, f)
    yi = ones(length(x)) .* y';
    yu = ones(length(x)) .* vcat(y[2:end],y[end])';
    yl = ones(length(y))' .* f;
    yr = ones(length(y))' .* vcat(f[2:end],f[end]);
    R = area_ratio.(yi,yu,yl,yr);
    return R
end

@doc raw"""
    area_ratio_logbins(kkh, kkz, logλz, f)

Same as [`area_ratio_grid`](@ref "area_ratio_grid") but for logarithmic grids. 

* `kkh`: horizontal wave vectors grid points
* `kkz`: vertical wave vectors grid points
* `logλz`: log(λz) where λz is the logarithmic increment of the vertical wave vectors grid
* `f`: values of `f` vs `kkh`

Note: convolution allows to integrate over `kkz` < `f`.

"""
function area_ratio_logbins(kkh, kkz, logλz, f)
    iff = zeros(length(kkh))
    @. iff = (log(f) - log(kkz[1])) / logλz + 1 
    yi = ones(length(kkh)) .* vcat(1:length(kkz))';
    yu = ones(length(kkh)) .* vcat(vcat(2:length(kkz)),[length(kkz)])';
    yl = ones(length(kkz))' .* iff;
    yr = ones(length(kkz))' .* vcat(iff[2:end],iff[end]);
    R = area_ratio.(yi,yu,yl,yr);

    return R
end