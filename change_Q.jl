"""Increasing Q (upper energy flux) leads to unimodal Pr[sustainability] response"""
"""Testing this observation in a series of σs"""

# work on this script towards Figure 2 in the paper. Will need to change the parameters to fit realistic ecosystems
# for example, d=100 instead of 1

include("../src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra
using Plots, ProgressMeter, IterTools

## --- emergence ---
# ec_1 = ecosys_config(K=8, S_type=:indiv, conne=1.0, k_param=0.2, d_param=:lognormal, seed=123);
# ec_1 = ecosys_config(K=4, S_type=:indiv, conne=1.0, k_param=0.2, d_param=:allometric, N0_param=:lognormal, seed=123);

# σs = generate_sigma_arrays(ec_1, 20);
# Q_range = Vector(1e-4:0.1:15.0);

# vols_all_1 = []; devols_all_1 = [];
# Random.seed!(345);
# for i in eachindex(σs)
#     @info "Compute Vol for $i in $(length(σs))"
#     p_i = generate_problem(ec_1, σs[i])
#     vols = volume_range_EFD(p_i, Q_range)
#     devols = [S^p_i.K / (factorial(p_i.K) * prod(p_i.N⁰)) for S in Q_range];
#     push!(vols_all_1, vols)
#     push!(devols_all_1, devols)
# end

## ---- persistence ---- 
ec = ecosys_config(K=4, S_type=:total, conne=1.0, k_param=0.1, d_param=:lognormal, seed=123);
σs = generate_sigma_arrays(ec, 20);
Q_range = Vector(1e-4:0.35:15.0);

vols_all = []; devols_all = [];
Random.seed!(345);
for i in eachindex(σs)
    @info "Compute Vol for $i in $(length(σs))"
    p_i = generate_problem(ec, σs[i])
    vols = volume_range_EFD(p_i, Q_range, n_sample=10^4)
    devols = [S^p_i.K / (factorial(p_i.K) * prod(p_i.N⁰)) for S in Q_range];
    push!(vols_all, vols)
    push!(devols_all, devols)
end

# Plotting
prob_pers = [vols./devols for (vols, devols) in zip(vols_all, devols_all)];
prob_pers = [replace(x -> isinf(x) ? 0.0 : x, curve) for curve in prob_pers];

peaks = [extract_peak(Q_range, curve) for curve in prob_pers];
P_opt = getindex.(peaks, 1); p_per = getindex.(peaks, 2);

plot(Q_range, prob_pers, xlabel="Total Energy Supply", ylabel="P[Persistence]", legend=false, color = :olivedrab, ylim=(0,0.5), linewidth=1.5, framestyle=:box, foreground_color_border=:black, borderlinewidth=1, grid=false)
scatter!(P_opt, p_per, label="Peak", color=:black, markersize=4, marker=:circle, legend=false)