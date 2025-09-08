
using EnerFeas
using Random, Distributions, LinearAlgebra
using Plots, ProgressMeter, IterTools
using LaTeXStrings

labels = ["σ₀=0.1", "σ₀=0.25", "σ₀=0.5", "σ₀=0.75", "σ₀=1.0"];
colors = [:green, :blue, :black, :red, :darkorange];

ecs = [ecosys_config(4, σ0=s, seed=42) for s in [0.1, 0.25, 0.5, 0.75, 1.0]];
vols_total = Vector{Vector{Float64}}(undef, 5); 

Random.seed!(345);
Q_range_t = vcat(1e-4:0.25:100.0);
for (i, ec) in enumerate(ecs_total)
    @info "compute the $(i) config"
    σ = generate_sigma_arrays(ec, 1);
    p = generate_model_system(ec, σ)
    vols_total[i] = volume_range_EFD(p, :matr, Q_range_t, n_sample=2*10^4)
end

devols_t = [Q^4/factorial(4) for Q in Q_range_t];

plt_t = plot(
    framestyle=:box, grid=false,
    xlabel="Upper Energy Bound " * L"(Q)",
    ylabel="Prob. Feasibility " * L"(\mathbb{P}^F)",
    xaxis = :log10,
    xlim = (1.0, 100.0),
    ylim = (0.0, 0.25),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    size=(225,180),
    legend=(0.15,0.93),
    foreground_color_legend = nothing,
    background_color = :transparent
);

for (v, la, co) in zip(vols_total, labels, colors)
    plot!(plt_t, Q_range_t, v ./ devols_t; color=co, linewidth=1.0, label=la)
end
display(plt_t)
savefig(plt_t, "figures/sigma0_total.pdf");

# ----- Individual Energy Bound -----

vols_indiv = Vector{Vector{Float64}}(undef, 5); 

Random.seed!(345);
Q_range_i = vcat(1e-4:0.25:1000.0);
for (i, ec) in enumerate(ecs)
    @info "compute the $(i) config"
    σ = generate_sigma_arrays(ec, 1);
    p = generate_model_system(ec, σ)
    vols_indiv[i] = volume_range_EFD(p, :init, Q_range_i)
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
    plot!(plt_i, Q_range_i, v ./ devols_i; color=co, linewidth=1.0, label=la)
end
display(plt_i)
savefig(plt_i, "figures/sigma0_indiv.pdf");

using JLD2
@save "data/sensitivity/sigma0.jld" ecs_total ecs_indiv Q_range_t Q_range_i vols_total vols_indiv