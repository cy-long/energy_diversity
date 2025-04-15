"Increasing P for persistence constraint generate unimodal response"

include("../src/EnerFeas.jl");
using .EnerFeas

using Random, Distributions, LinearAlgebra
using Plots, StatsPlots, DataFrames, ProgressMeter, IterTools

# --- emergence ---
# ec_1 = ecosys_config(K=8, S_type=:indiv, conne=1.0, k_param=0.2, d_param=:lognormal, seed=123);
ec_1 = ecosys_config(K=4, S_type=:indiv, conne=1.0, k_param=0.2, d_param=:allometric, N0_param=:lognormal, seed=123);

σs = generate_sigma_arrays(ec_1, 20);
P_range = Vector(1e-4:0.1:15.0);

vols_all_1 = []; devols_all_1 = [];
Random.seed!(345);
for i in eachindex(σs)
    @info "Compute Vol for $i in $(length(σs))"
    p_i = generate_problem(ec_1, σs[i])
    vols = volume_range_EFD(p_i, P_range)
    devols = [S^p_i.K / (factorial(p_i.K) * prod(p_i.N⁰)) for S in P_range];
    push!(vols_all_1, vols)
    push!(devols_all_1, devols)
end


# ---- persistence ----
# ec = ecosys_config(K=4, S_type=:total, conne=1.0, k_param=0.2, d_param=:lognormal, seed=123);
ec = ecosys_config(K=4, S_type=:total, conne=1.0, k_param=0.05, d_param=:allometric, N0_param=:lognormal, seed=123);

vols_all = []; devols_all = [];
Random.seed!(345);
for i in eachindex(σs)
    @info "Compute Vol for $i in $(length(σs))"
    p_i = generate_problem(ec, σs[i])
    vols = volume_range_EFD(p_i, P_range)
    devols = [S^p_i.K / (factorial(p_i.K) * prod(p_i.N⁰)) for S in P_range];
    push!(vols_all, vols)
    push!(devols_all, devols)
end

# Plotting
prob_pers = [vols./devols for (vols, devols) in zip(vols_all_1, devols_all_1)];
plot(P_range, prob_pers, xlabel="Total Energy Supply", ylabel="P[Emergence]", legend=false, color = :olivedrab, ylim=(0,1.0), linewidth=1.5, framestyle=:box, foreground_color_border=:black, borderlinewidth=2, grid=false)

prob_pers = [vols./devols for (vols, devols) in zip(vols_all, devols_all)];
plot(P_range, prob_pers, xlabel="Total Energy Supply", ylabel="P[Persistence]", legend=false, color = :olivedrab, ylim=(0,0.1), linewidth=1.5, framestyle=:box, foreground_color_border=:black, borderlinewidth=2, grid=false)