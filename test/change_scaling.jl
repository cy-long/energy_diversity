"""Instrinsic scaling in all energy related variables"""

include("../src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra
using Plots, ProgressMeter, IterTools
using LaTeXStrings

ecs_total = [ecosys_config(K=4,S_type=:total,k_param=0.0,d_param=s,n_scale=s,seed=42) for s in [1.0, 10.0, 100.0]];

Q_range_l = vcat(1.0:0.2:10.0, 10.0:1.0:100.0, 100.0:5.0:2000.0);
Q_range_h = vcat(1.0, 100.0:10.0:2000.0, 2000.0:50.0:10000.0);
Q_ranges_t = [Q_range_l, Q_range_l, Q_range_h];

vols_total = Vector{Vector{Vector{Float64}}}(undef, 3);
Random.seed!(345);

for (i, ec) in enumerate(ecs_total)
    σs = generate_sigma_arrays(ec, 5);
    vols = Vector{Vector{Float64}}(undef, 5)
    @info "compute the $(i) config"
    for (j, σ) in enumerate(σs)
        p = generate_problem(ec, σ)
        vols[j] = volume_range_EFD(p, Q_ranges_t[i], n_sample=2*10^4)  # this is a Vector
    end
    vols_total[i] = vols
end

plt_t = plot(
    framestyle=:box, grid=false,
    xlabel="Upper Energy Bound " * L"(Q)",
    ylabel="Prob. Feasibility " * L"(\mathbb{P}^F)",
    xaxis = :log10,
    xlim = (1.0,10000.0),
    ylim = (0.0,0.2),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    size=(225,180),
    legend=(0.85,0.93),
    foreground_color_legend = nothing,
    legend_border = false,
    background_color = :transparent,
    lw=1.0,
    );

labels = ["k=10⁰","k=10¹","k=10²"]; colors = [:green, :blue, :red];

for (vols, label, color, Q_range) in zip(vols_total, labels, colors, Q_ranges_t)
    devols = [Q^4/factorial(4) for Q in Q_range];
    for (i, v) in enumerate(vols)
        plot!(
            plt_t, Q_range, v ./ devols;
            color=color, linewidth=1.0,
            label=(i == 1 ? label : "")
        )
    end
end

display(plt_t)
savefig(plt_t, "figures/scaling_total.pdf");


# ----- Individual Energy Bound -----
ecs_indiv = [ecosys_config(K=4,S_type=:indiv,k_param=0.0,d_param=s,n_scale=s,seed=42) for s in [1.0, 10.0, 100.0]];
Q_ranges_i = [vcat(1.0:1.0:10000.0) for _ in 1:3];
vols_indiv = Vector{Vector{Vector{Float64}}}(undef, 3);
Random.seed!(345);

for (i, ec) in enumerate(ecs_indiv)
    σs = generate_sigma_arrays(ec, 5);
    vols = Vector{Vector{Float64}}(undef, 5)
    @info "compute the $(i) config"
    for (j, σ) in enumerate(σs)
        p = generate_problem(ec, σ)
        vols[j] = volume_range_EFD(p, Q_ranges_i[i], n_sample=2*10^4)  # this is a Vector
    end
    vols_indiv[i] = vols
end

plt_i = plot(
    framestyle=:box, grid=false,
    xlabel="Upper Energy Bound " * L"(Q)",
    ylabel="Prob. Viability " * L"(\mathbb{P}^V)",
    xaxis = :log10,
    xlim = (1.0,10000.0),
    ylim = (0.0,0.75),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    size=(225,180),
    legend=(0.1,0.93),
    foreground_color_legend = nothing,
    legend_border = false,
    background_color = :transparent,
    lw=1.0,
    );

labels = ["k=10⁰","k=10¹","k=10²"]; colors = [:green, :blue, :red];
for (vols, label, color, Q_range) in zip(vols_indiv, labels, colors, Q_ranges_i)
    devols = [Q^4/factorial(4) for Q in Q_range];
    for (i, v) in enumerate(vols)
        plot!(
            plt_i, Q_range, v ./ devols;
            color=color, linewidth=1.0,
            label=(i == 1 ? label : "")
        )
    end
end

display(plt_i)
savefig(plt_i, "figures/scaling_indiv.pdf");