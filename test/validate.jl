"""Verify the volume calculation, test for best hyperparameters (such as n_sample), with changing K"""

include("../src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra
using ProgressMeter, Plots, SpecialFunctions

# testing volume of the intersection of ∑ᵢxᵢ²≤r² and the positive orthant of xᵢ≥0
Random.seed!(1234)
function vol_spherical_simplex(r::Float64, K::Int)
    return (pi^(K/2) * r^K) / ((gamma(K/2 + 1) * 2^K))
end

K_range = [2,4,8,12];
vols_stand = Dict{Int, Vector{Float64}}();
vols_samp = Dict{Int, Vector{Float64}}();

for K in K_range
    @info "K = $K"
    r_range = exp.(range(0, log(100), length=25));
    vol_stan_K = [vol_spherical_simplex(r, K) for r in r_range];

    A = -Matrix{Float64}(I, (K, K));
    b = zeros(K);
    regions = [InterPolyBalls(A, b, [Sphere(zeros(K), r)], :total) for r in r_range];

    vol_samp_K = zeros(Float64, length(regions))
    @showprogress for (i, reg) in pairs(regions)
        vol_samp_K[i] = volume_domain(reg, 10, 2, false)
    end

    vols_stand[K] = vol_stan_K;
    vols_samp[K] = vol_samp_K;
end

# plot the results: radius and size vs. relative error
plt = plot(
    xaxis=:log10, yaxis=:log10, 
    xlabel="Radius r", 
    ylabel="Relative Error", 
    legend=:bottomright,
    grid=false,
    yticks=[1e-4, 1e-3, 1e-2, 1e-1],size=(400,300));

for K in K_range
    r_range = exp.(range(0, log(100), length=25));
    vol_stan_K = vols_stand[K];
    vol_samp_K = vols_samp[K];

    rel_error = abs.(vol_samp_K .- vol_stan_K) ./ vol_stan_K;
    @info "$K has mean_rel_error: $(mean(rel_error))"
    plot!(plt, r_range, rel_error, label="K=$K", marker=:circle, markersize=3, markerstrokecolor=:auto)
end

savefig(plt, "figures/verify_volume.svg")