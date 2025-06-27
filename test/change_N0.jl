"""Changing the average level of N⁰"""

include("../src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra
using Plots, ProgressMeter, IterTools
using LaTeXStrings

labels = ["N⁰=1", "N⁰=√5", "N⁰=5"]; colors = [:green, :blue, :red]; N0s = [1.0, sqrt(5), 5.0];

# ----- Total Energy Bound -----

ecs_total = [ecosys_config(K=4, S_type=:total, N0_param = N0, k_param=0.1, seed=42) for N0 in [1.0, sqrt(5), 5.0]];

Q_range_1 = vcat(1.0:0.2:50.8, 52.0:3.8:1000.0);
Q_range_2 = vcat(1.0:0.75:200.0, 200.0:3.42:1000.0);
Q_range_3 = vcat(1.0:1.0:100.0, 102.25:2.25:1000.0);
Q_ranges_t = [Q_range_1, Q_range_2, Q_range_3];

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
    xlim = (1.0,1000.0),
    ylim = (0.0, 0.3),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    size=(225,180),
    legend=(0.15,0.93),
    foreground_color_legend = nothing,
    background_color = :transparent
    );

for (vols, la, co, N0, Q_range) in zip(vols_total, labels, colors, N0s, Q_ranges_t)
    devols = [(Q/N0)^4/factorial(4) for Q in Q_range];
    for (i, v) in enumerate(vols)
        plot!(
            plt_t, Q_range, v ./ devols;
            color=co, linewidth=1.0,
            label=(i == 1 ? la : "")
        )
    end
end

display(plt_t)
savefig(plt_t, "figures/N0_total.pdf")

# ----- Individual Energy Bound -----

ecs_indiv = [ecosys_config(K=4, S_type=:indiv, N0_param = N0, k_param=0.1, seed=42) for N0 in [1.0, sqrt(5), 5.0]];
Q_ranges_i = [vcat(1.0:0.5:1000.0) for _ in 1:3];
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
    xlim = (1.0,1000.0),
    ylim = (0.0, 0.65),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    size=(225,180),
    legend=(0.15,0.93),
    foreground_color_legend = nothing,
    background_color = :transparent,
    );

for (vols, label, color, N0, Q_range) in zip(vols_indiv, labels, colors, N0s, Q_ranges_i)
    devols = [(Q/N0)^4/factorial(4) for Q in Q_range];
    for (i, v) in enumerate(vols)
        plot!(
            plt_i, Q_range, v ./ devols;
            color=color, linewidth=1.0,
            label=(i == 1 ? label : "")
        )
    end
end

display(plt_i)
savefig(plt_i, "figures/N0_indiv.pdf")