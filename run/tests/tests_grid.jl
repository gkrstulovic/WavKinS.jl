# Tests for grid

push!(LOAD_PATH, "../")
using WavKinS
using GLMakie

println("---------------------------------------------------------------------")
println("Test area_ratio and area_ratio_logbins functions")
println("")

# Check area_ratio (with lin bins)
for M=3:10
    x = LinRange(0,2*pi,2^M);
    y = LinRange(0,2,2^(M+1));
    dx = x[2]-x[1];
    dy = y[2]-y[1];
    
    #f = pi*ones(1,length(x))/2;
    #Itheo=pi^2;
    
    #f = 0.1 * x;
    #Itheo = 0.1 * (2*pi)^2 / 2;

    #f = 0.02 * x.^2;
    #Itheo = 0.02 * (2*pi)^3 / 3;

    f = ones(length(x)) + sin.(x);
    Itheo=2*pi;
    
    R = area_ratio_grid(x,y,f);
    I = sum(R, dims=[1,2]) * dx * dy;
    I = I[1,1]
    error_inte = abs((I-Itheo)/Itheo);
    println(2^M, "   ", error_inte)
    fig, ax = plot_heatmap_base!(x, y, R, title="Test area_ratio_grid", xscale=identity, yscale=identity, zscale=identity)
    display(fig)
    sleep(5)
end


# Check area_ratio_logbins
for M=3:10
    Mh = 2^M
    Mz = 2^M
    khmin = 1e-4
    khmax = 2*pi
    kzmin = 1e-4
    kzmax = 2.0
    Nk = wave_spectrum_khkz(khmin,khmax,Mh,kzmin,kzmax,Mz)

    kkh = Nk.kkh
    kkz = Nk.kkz
    logλz = Nk.logλz
    
    #f = pi*ones(1,length(kkh))/2;
    #Itheo=pi^2;
    
    f = 0.1 * kkh;
    Itheo = 0.1 * (2*pi)^2 / 2;

    #f = 0.02 * kkh.^2;
    #Itheo = 0.02 * (2*pi)^3 / 3;

    #f = ones(Mh) + sin.(kkh);
    #Itheo=2*pi;
    
    R = area_ratio_logbins(kkh,kkz,logλz,f);
    @. Nk.nk = R
    integ = WavKinS.integrate_with_log_bins_khkz()
    I = integrate(integ, Nk, 1, -1, 1, -1)
    error_inte = abs((I-Itheo)/Itheo);
    println(Mh, "   ", error_inte)
    fig, ax = plot_heatmap_base!(kkh, kkz, R, title="Test area_ratio_logbins", xscale=identity, yscale=identity, zscale=identity);
    display(fig)
    sleep(5)
end
