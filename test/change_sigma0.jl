include("../src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra
using Plots, ProgressMeter, IterTools
using LaTeXStrings

labels = ["σ₀=0.1", "σ₀=0.25", "σ₀=0.75", "σ₀=1.0"]; colors = [:green, :blue, :red, :darkorange];

ecs_total = [ecosys_config(K=4, S_type=:total, k_param=0.1, seed=42, ϵ=a) for a in [0.1, 0.25, 0.75, 1.0]];
vols_total = Vector{Vector{Float64}}(undef, 4); 

Random.seed!(345);
Q_range_t = vcat(1e-4:0.25:100.0);
for (i, ec) in enumerate(ecs_total)
    @info "compute the $(i) config"
    σ = generate_sigma_arrays(ec, 1);
    p = generate_problem(ec, σ)
    vols_total[i] = volume_range_EFD(p, Q_range_t, n_sample=2*10^4)
end

devols_t = [Q^4/factorial(4) for Q in Q_range_t];

plt_t = plot(
    framestyle=:box, grid=false,
    xlabel="Upper Energy Bound " * L"(Q)",
    ylabel="Prob. Feasibility " * L"(\mathbb{P}^F)",
    xaxis = :log10,
    xlim = (1.0,100.0),
    # ylim = (0.0, 0.3),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    size=(225,180),
    legend=(0.15,0.93),
    foreground_color_legend = nothing,
    background_color = :transparent
);

for (v, la, co) in zip(vols_total, labels, colors)
    plot!(plt_t, Q_range_t, v ./ devols_t; color=co, linewidth=1.5, label=la)
end
display(plt_t)
savefig(plt_t, "figures/sigma0_total.pdf");


ecs_indiv = [ecosys_config(K=4, S_type=:indiv, k_param=0.1, seed=42, ϵ=a) for a in [0.1, 0.25, 0.75, 1.0]];
vols_indiv = Vector{Vector{Float64}}(undef, 4); 

Random.seed!(345);
Q_range_i = vcat(1e-4:0.25:1000.0);
for (i, ec) in enumerate(ecs_indiv)
    @info "compute the $(i) config"
    σ = generate_sigma_arrays(ec, 1);
    p = generate_problem(ec, σ)
    vols_indiv[i] = volume_range_EFD(p, Q_range_i)
end


devols_i = [Q^4/factorial(4) for Q in Q_range_i];


plt_i = plot(
    framestyle=:box, grid=false,
    xlabel="Upper Energy Bound " * L"(Q)",
    ylabel="Prob. Viability " * L"(\mathbb{P}^V)",
    xaxis = :log10,
    xlim = (1.0,100.0),
    # ylim = (0.0, 0.3),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    size=(225,180),
    legend=(0.15,0.93),
    foreground_color_legend = nothing,
    background_color = :transparent
);

for (v, la, co) in zip(vols_indiv, labels, colors)
    plot!(plt_i, Q_range_i, v ./ devols_i; color=co, linewidth=1.5, label=la)
end
display(plt_i)
savefig(plt_i, "figures/sigma0_indiv.pdf");

using JLD2
@save "data/change_eps.jld2" ecs_total ecs_indiv Q_range_t Q_range_i vols_total vols_indiv

# Let's then understand the relationship between ϵ and self-regulation
# using LinearAlgebra
# using Plots

N = 20000; K = 4; n_var = 1.0;
As = [rand(Normal(0, n_var),K, K) for _ in 1:N];
Ps = [(a + a') / 2 for a in As];
cs = [-minimum(eigvals(p), init=0.0) for p in Ps];
# historgram of cs
using StatsBase
histogram(cs, bins=100, xlabel="Minimum Eigenvalue", ylabel="Frequency", normalize=true)


# plt = plot(aspect_ratio = :equal)
# ϵs = [0.05, 0.1, 0.25, 0.5, 0.75, 1.0, 1.25];

# for ϵ in ϵs
#     Random.seed!(10);
#     σ = rand(Normal(0,3.0), 2, 2);
#     # c = minimum(eigvals((σ + σ') / 2), init=0.0)
#     σ[diagind(σ)] .+= ϵ;

#     Λ = inv(σ);
#     P = (Λ + Λ') / 2

#     # visualize the ellipsoid
#     L = Matrix(cholesky(P).U);
#     pts = [[cos.(θ),sin.(θ)] for θ in 0:0.05:2π];
#     pts = [inv(L) * pt for pt in pts];

#     plot!(plt, [p[1] for p in pts], [p[2] for p in pts], alpha=0.35, linewidth=3, label=@sprintf("ϵ=%.2f", ϵ))
# end

# display(plt)

# rms = [randn(10,10) for _ in 1:100000]
# eigs = [eigvals((a+a')/2) for a in rms]
# eigmins = [minimum(e) for e in eigs]
# histogram(eigmins, bins=100, xlabel="Eigenvalues", ylabel="Frequency", title="Min Eigvals of random symmetric mats")