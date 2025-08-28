""" Changing the energy exchange network of the system"""

include("../src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra
using Plots, ProgressMeter, IterTools
using LaTeXStrings

labels = [L"s_\sigma=0.5", L"s_\sigma=1.0", L"s_\sigma=2.0"]; colors = [:green, :blue, :red]; scs = [0.5, 1.0, 2.0];

ecs_total = [ecosys_config(K=4, S_type=:total, n_scale=sc, ϵ_param=0.0, seed=42) for sc in scs];

Q_range_t = vcat(1e-4:0.2:100.0);

vols_total = Vector{Vector{Vector{Float64}}}(undef, 3); 
Random.seed!(345);

for (i, ec) in enumerate(ecs_total)
    σs = generate_sigma_arrays(ec, 1);
    vols = Vector{Vector{Float64}}(undef, 1)
    @info "compute the $(i) config"
    # for (j, σ) in enumerate(σs)
    j = 1; σ = σs
    p = generate_problem(ec, σ)
    vols[j] = volume_range_EFD(p, Q_range_t, n_sample=2*10^4)  # this is a Vector
    # end
    vols_total[i] = vols
end

plt_t = plot(
    framestyle=:box, grid=false,
    xlabel="Upper Energy Bound " * L"(Q)",
    ylabel="Prob. Realization " * L"(\mathbb{P}^R)",
    xaxis = :log10,
    xlim = (1.0, 100.0),
    # ylim = (0.0, 0.35),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    # size=(225,180),
    size = (400,300),
    legend=(0.8,0.93),
    foreground_color_legend = nothing,
    background_color = :transparent
);

for (vols, la, co, sc) in zip(vols_total, labels, colors, scs)
    devols = [(Q/1.0)^4/factorial(4) for Q in Q_range_t];
    for (i, v) in enumerate(vols)
        plot!(
            plt_t, Q_range_t, v ./ devols;
            color=co, linewidth=1.0,
            label=(i == 1 ? la : "")
        )
    end
end

display(plt_t)
# savefig(plt_t, "figures/sigma_total.pdf");

# ----- Individual ------
ecs_indiv = [ecosys_config(K=4, S_type=:indiv, n_scale=sc, ϵ_param=0.0, seed=42) for sc in scs];

Q_range_i = vcat(1e-4:0.2:1000.0);

vols_indiv = Vector{Vector{Vector{Float64}}}(undef, 3); 
Random.seed!(345);

for (i, ec) in enumerate(ecs_indiv)
    σs = generate_sigma_arrays(ec, 5);
    vols = Vector{Vector{Float64}}(undef, 5)
    @info "compute the $(i) config"
    for (j, σ) in enumerate(σs)
        p = generate_problem(ec, σ)
        vols[j] = volume_range_EFD(p, Q_range_i, n_sample=2*10^4)  # this is a Vector
    end
    vols_indiv[i] = vols
end

plt_i = plot(
    framestyle=:box, grid=false,
    xlabel="Upper Energy Bound " * L"(Q)",
    ylabel="Prob. Initiation " * L"(\mathbb{P}^I)",
    xaxis = :log10,
    xlim = (1.0, 200.0),
    # ylim = (0.0, 0.35),
    guidefont=font(6),
    tickfont=font(5),
    legendfont=font(5),
    size=(225,180),
    # legend=(0.13,0.93),
    legend=false,
    foreground_color_legend = nothing,
    background_color = :transparent
);

for (k, vols, la, co, N0) in zip([0,1,2], vols_indiv, labels, colors, d0s)
    devols = [(Q/1.0)^4/factorial(4) for Q in Q_range_i];
    for (i, v) in enumerate(vols)
        plot!(
            plt_i, Q_range_i, v ./ devols .+ k*0.006;
            color=co, linewidth=1.0,
            label=(i == 1 ? la : "")
        )
    end
end

display(plt_i)
savefig(plt_i, "figures/sigma_indiv.pdf");