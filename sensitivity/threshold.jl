

using EnerFeas
using Random, Distributions, LinearAlgebra
using Plots, ProgressMeter, IterTools
using LaTeXStrings

labels = ["ϵ=0","ϵ=10⁻³","ϵ=10⁻¹"]; colors = [:green, :blue, :red]; thres = [0.0, 0.001, 0.1];

# ----- Total Energy Bound -----
ecs = [ecosys_config(4, ϵ=e, seed=42) for e in thres];
Q_range_t = vcat(1e-4:0.1:10, 10.25:0.3:100.0);

vols_total = Vector{Vector{Vector{Float64}}}(undef, 3); 
Random.seed!(345);

for (i, ec) in enumerate(ecs)
    if i > 2 continue end  # skip the first one, which is the baseline
    σs = generate_sigma_arrays(ec, 5);
    vols = Vector{Vector{Float64}}(undef, 5)
    @info "compute the $(i) config"
    for (j, σ) in enumerate(σs)
        p = generate_model_system(ec, σ)
        vols[j] = volume_range_EFD(p, :matr, Q_range_t, n_sample=2*10^4)
    end
    vols_total[i] = vols
end

plt_t = plot(
    framestyle=:box, grid=false,
    xlabel="Upper Energy Bound " * L"(Q)",
    ylabel="Prob. Feasibility " * L"(\mathbb{P}^F)",
    xaxis = :log10,
    xlim = (1.0,100.0),
    ylim = (0.0, 0.2),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    size=(225,180),
    legend=(0.82,0.93),
    foreground_color_legend = nothing,
    background_color = :transparent
);

lines = [:solid, :dot, :dot]; widths = [1.5, 1.0, 1.0];
for (vols, la, co, li, wi) in zip(vols_total, labels, colors, lines, widths)
    devols = [Q^4/factorial(4) for Q in Q_range_t];
    for (i, v) in enumerate(vols)
        plot!(
            plt_t, Q_range_t, v ./ devols;
            color=co, linewidth=wi, linestyle=li,
            label=(i == 1 ? la : "")
        )
    end
end

display(plt_t)
savefig(plt_t, "figures/thre_total.pdf");


# ----- Individual Energy Bound -----
Q_range_i = vcat(1e-4:0.1:100.0, 100.0:1.0:1000.0);
vols_indiv = Vector{Vector{Vector{Float64}}}(undef, 3);
Random.seed!(345);

for (i, ec) in enumerate(ecs)
    σs = generate_sigma_arrays(ec, 5);
    vols = Vector{Vector{Float64}}(undef, 5)
    @info "compute the $(i) config"
    for (j, σ) in enumerate(σs)
        p = generate_model_system(ec, σ)
        vols[j] = volume_range_EFD(p, :init, Q_range_i, n_sample=2*10^4)
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
    background_color = :transparent
);

lines = [:solid, :dot, :dot]; widths = [1.5, 1.0, 1.0];
for (vols, la, co, li, wi) in zip(vols_indiv, labels, colors, lines, widths)
    devols = [Q^4/factorial(4) for Q in Q_range_i];
    for (i, v) in enumerate(vols)
        plot!(
            plt_i, Q_range_i, v ./ devols;
            color=co, linewidth=wi, linestyle=li,
            label=(i == 1 ? la : "")
        )
    end
end

display(plt_i)
savefig(plt_i, "figures/thre_indiv.pdf");