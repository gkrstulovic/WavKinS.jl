# Tests for misc

push!(LOAD_PATH, "../")
using WavKinS
using GLMakie

println("---------------------------------------------------------------------")
println("Test cardano function")
println("")

Ip = 0:1000
Iq = 0:1000
pp = 10.0 .* Ip
qq = -1e-3 .* Iq
for p in pp, q in qq
    x = WavKinS.cardano(p,q)
    f = x^3 + p*x + q
    #if f > 1e-10
    #    println(f)
    #end
    @assert f < 1e-10
end
println("cardano function ok")


println("")
println("---------------------------------------------------------------------")


println("---------------------------------------------------------------------")
println("Test NewtonRaphson function")
println("")

f(x) = x^x - 100
fp(x) = x^x * (log(x) + 1)

x0 = NewtonRaphson(f, 1.0)
println("NewtonRaphson: ", x0, "  ", f(x0))
x0bis = NewtonRaphson_with_derivative(f, fp, 1.0)
println("NewtonRaphson_with_derivative: ", x0bis, "  ", f(x0bis))
x = vcat(0:0.001:5)

fig, ax = plot_1D_base!(x, f.(x), title="Test NewtonRaphson", xscale=identity, yscale=identity)
scatter!(ax,[x0],[0], label="num")
scatter!(ax,[x0bis],[0], label="with derivative")
xlims!(3,4)
ylims!(-1,1)
axislegend(ax, position = :lt)
display(fig)
sleep(5)

# Check for Acoustic2D_khkz

ω = WavKinS.ω_Acoustic_disp
k = 0.365
k1 = LogRange(1e-2,1e1,256)
k2 = LogRange(1e-2,1e1,256)

K1, K2 = meshgrid(k1, k2)
K = k*ones(size(K1))
dω = ω.(K) - ω.(K1) - ω.(K2)

k2sol = zeros(length(k1))

for i in eachindex(k2sol)
    f(x) = ω(k) - ω(k1[i]) - ω(x)
    k2sol[i] = NewtonRaphson(f, k1[i])
end

#fig, ax = plot_heatmap_base!(k1,k2,dω; zscale=identity, zlims=(-10,10), colormap=:coolwarm)
fig, ax = plot_surface_base!(k1,k2,dω; xlabel=L"$k_1$", ylabel=L"$k_2$", zlabel=L"$\omega_k - \omega_1 - \omega_2$", title="NewtonRaphson for Acoustic", xscale=log10, zscale=identity, zlims=(-1,1))
lines!(ax, log10.(k1[k2sol.>0]), log10.(k2sol[k2sol.>0]), color="magenta", lw=2, legend=false)
display(fig)

println("")
println("---------------------------------------------------------------------")

