"""Explore the relationship between σ and P-p curves"""

include("../src/EnerFeas.jl")
using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots
using Loess
using Printf

P_range = Vector(1e-4:0.25:20.0);
scales = [0.01, 0.1, 1.0];
# vars = [0.5, 1.0, 2.0];
colors = [:red, :green, :blue];

Random.seed!(345);
plt = plot(framestyle=:box, grid=false, xlim=(0,20), ylim=(0,0.6));
for i in eachindex(scales)
    ec = ecosys_config(K=4, S_type=:total, n_scale=scales[i], n_var=1.0, conne=1.0, k_param=0.2, d_param=0.1, N0_param=1.0, seed=50)
    σ = generate_sigma_arrays(ec, 10);
    all_vols = []; all_devols = [];
    for j in [1]
        @info "Compute Volume in $j"
        p = generate_problem(ec, σ[j]);
        vols = volume_range_EFD(p, P_range);
        devols = [S^p.K / (factorial(p.K) * prod(p.N⁰)) for S in P_range];
        push!(all_vols, vols)
        push!(all_devols, devols)
    end
    prob_pers = [vols./devols for (vols, devols) in zip(all_vols, all_devols)];
    prob_pers = [replace(x -> isinf(x) ? 0.0 : x, curve) for curve in prob_pers];
    plot!(plt, P_range, prob_pers, color = colors[i], alpha=0.2, label=false, linewidt=1.5)
    peaks = [extract_peak(P_range, curve) for curve in prob_pers];
    P_opt = getindex.(peaks, 1); p_per = getindex.(peaks, 2);
    scatter!(plt, P_opt, p_per, label=@sprintf("scales=%.2f", scales[i]), xlabel="P_opt", ylabel="max Pr", legend=:topright, color = colors[i])
end
display(plt)
