"""Define and measure the score of saturation / unimodal response"""

include("../src/EnerFeas.jl")
using .EnerFeas, JLD2, Distributions, Random, Plots, StatsPlots, DataFrames

##1 demonstrate the unimodal score
ec = ecosys_config(K=4, S_type=:total, conne=1.0, ϵ_param=0.2, d_param=0.1, N0_param=1.0, seed=100);
σ = generate_sigma_arrays(ec, 50);
Q_range = Vector(1e-4:0.1:20.0);

all_vols = []; all_devols = []; scores = [];

Random.seed!(345);
for i in eachindex(σ)
    @info "Compute Volume in $i"
    p = generate_problem(ec, σ[i]);
    vols = volume_range_EFD(p, Q_range,n_sample=10^4);
    devols = [S^p.K / (factorial(p.K) * prod(p.N⁰)) for S in Q_range];
    push!(all_vols, vols)
    push!(all_devols, devols)
    push!(scores, score_unimodal(Q_range, vols ./ devols))
end

# let's generate a unbiased randomwalk in [0,20]
scores_null = [];
for i in 1:50
    dy = rand(Normal(0,1),length(Q_range)-1);
    xx = Q_range[2:end];
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


##2 correlation of characteristics of peak location and magnitude
peaks = [extract_peak(Q_range, vols ./ devols) for (vols, devols) in zip(all_vols, all_devols)];
P_opt = getindex.(peaks, 1); p_per = getindex.(peaks, 2);
scatter(P_opt, p_per, label="K=4", xlabel="P_opt", ylabel="max Pr", legend=:topright, xlim=(0, 20), ylims=(0, 0.5));

##3 correlation of characteristics in persistence and emergence
# generate data for linear problems
ec1 = ecosys_config(K=4, S_type=:indiv, conne=1.0, ϵ_param=0.2, d_param=0.1, N0_param=1.0, seed=100);
# generate additional individual problems
all_vols1 = []; all_devols1 = [];
for i in eachindex(σ)
    @info "Compute Volume in $i"
    p = generate_problem(ec1, σ[i]);
    vols = volume_range_EFD(p, Q_range);
    devols = [S^p.K / (factorial(p.K) * prod(p.N⁰)) for S in Q_range];
    push!(all_vols1, vols)
    push!(all_devols1, devols)
end

p_eme = [extract_peak(Q_range, vols ./ devols)[2] for (vols, devols) in zip(all_vols1, all_devols1)];
p_per = [extract_peak(Q_range, vols ./ devols)[2] for (vols, devols) in zip(all_vols, all_devols)];
p_per[isnan.(p_per)] .= 0.0;

scatter(p_eme, p_per, label="K=4",  xlabel="P_emergence", ylabel="P_persistence", legend=:topright, xlim=(0, 1), ylims=(0, 1))
cor(p_eme, p_per)