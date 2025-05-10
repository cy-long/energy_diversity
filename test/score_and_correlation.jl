"""Define and measure the score of saturation / unimodal response"""

include("../src/EnerFeas.jl")
using .EnerFeas, JLD2, Distributions, Random, Plots, StatsPlots, DataFrames

##1 demonstrate the unimodal score
ec = ecosys_config(K=4, S_type=:total, conne=1.0, k_param=0.2, d_param=0.1, N0_param=1.0, seed=100);
σ = generate_sigma_arrays(ec, 50);
P_range = Vector(1e-4:0.1:20.0);

all_vols = []; all_devols = []; scores = [];

Random.seed!(345);
for i in eachindex(σ)
    @info "Compute Volume in $i"
    p = generate_problem(ec, σ[i]);
    vols = volume_range_EFD(p, P_range,n_sample=10^4);
    devols = [S^p.K / (factorial(p.K) * prod(p.N⁰)) for S in P_range];
    push!(all_vols, vols)
    push!(all_devols, devols)
    push!(scores, score_unimodal(P_range, vols ./ devols))
end

# let's generate a unbiased randomwalk in [0,20]
scores_null = [];
for i in 1:50
    dy = rand(Normal(0,1),length(P_range)-1);
    xx = P_range[2:end];
    yy = cumsum(dy);
    push!(scores_null, score_unimodal(xx, yy))
end

grouped_scores = [scores[.!isnan.(scores)]; scores_null];
labels = [repeat(["Observed"], length(scores[.!isnan.(scores)])); repeat(["Random"], length(scores_null))];

df = DataFrame(score = grouped_scores, group = labels);
plt_score = plot(framestyle=:box, grid=false);
@df df boxplot!(plt_score, :group, :score,
    linewidth = 1.5,
    legend = false,
    outliercolor = :black,
    outliershape = :circle,
    title = "Unimodal Scores");
display(plt_score)