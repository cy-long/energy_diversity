include("../src/EnerFeas.jl");
using .EnerFeas, LinearAlgebra, Random, Plots, JLD2

colors = [
    RGB(0.99,0.55,0.38),  # salmon
    RGB(1.00,0.85,0.18),  # gold
    RGB(0.65,0.85,0.33),  # lime green
    RGB(0.40,0.76,0.65),  # teal
    RGB(0.55,0.63,0.80),  # soft blue
    RGB(0.91,0.54,0.76)   # lavender
];

K_range = [3,4,5,6,7,8];
P_range = Vector(0.01:1.0:300.0);

# case1 random σ
# Random.seed!(1010);
Random.seed!(1989);
σ = generate_sigma(K_range[end], false, true, 1.0, 2.0, 0.85)

# case2 energetically constrained σ
Random.seed!(999);
σ = generate_sigma(K_range[end], true, true, 1.0, 2.0, 0.85)

# case3 trophic chain
σ = generate_trophic(K_range[end], 1.0, 1.0, 0.5, 0.1)

## --- compare different K and their curves for emergence ---

vols_all = Vector{Vector{Float64}}(undef, length(K_range));
devols_all = Vector{Vector{Float64}}(undef, length(K_range));
for (i,K) in pairs(K_range)
    @info "In K = $K"
    ec_k = ecosys_config(K=K, S_type=:indiv, k_param=0.05, d_param=1.0, N0_param=1.5, seed=100)
    σ_k = σ[1:K, 1:K]
    if eigvals(σ_k + σ_k') |> minimum < 0.0
        throw(ArgumentError("submatrix σ_k is not pos. def."))
    end
    p_k = generate_problem(ec_k, σ_k)
    vols_all[i] =  volume_range_EFD(p_k, P_range)
    devols_all[i] = [S^p_k.K / (factorial(p_k.K) * prod(p_k.N⁰)) for S in P_range]
end

plt_em = plot(framestyle=:box, grid=false, legend=:topleft,
    xlabel="Max. Total Supply", ylabel="Prob. to sustain baseline biomass", fontsize=12, size=(450, 360));

for (i, K) in pairs(K_range)
    curve = vols_all[i] ./ devols_all[i]
    idx = findall(curve .> 0.0)
    scatter!(plt_em, [P_range[idx[1]]], [curve[idx[1]]], mode="markers", color=colors[i], marker=:square, label="", markersize=4, markerstrokecolor=:auto)
    plot!(plt_em, P_range[idx], curve[idx], label="K=$K", lw=1.5, xaxis=:log10, color=colors[i])
end
plot!(plt_em, title="Trophic σ")
savefig(plt_em, "figures/compare_K_emer_trop.pdf");
@save "data/compare_K_emer_trop.jld" K_range P_range σ vols_all devols_all plt_em

## --- compare different K and their curves for persistence ---
vols_all = Vector{Vector{Float64}}(undef, length(K_range));
devols_all = Vector{Vector{Float64}}(undef, length(K_range));
for (i,K) in pairs(K_range)
    @info "In K = $K"
    ec_k = ecosys_config(K=K, S_type=:total, k_param=0.05, d_param=1.0, N0_param=1.5, seed=100)
    σ_k = σ[1:K, 1:K]
    if eigvals(σ_k + σ_k') |> minimum < 0.0
        throw(ArgumentError("submatrix σ_k is not pos. def."))
    end
    p_k = generate_problem(ec_k, σ_k)
    vols_all[i] =  volume_range_EFD(p_k, P_range, n_sample=10^5, n_layer=20)
    devols_all[i] = [S^p_k.K / (factorial(p_k.K) * prod(p_k.N⁰)) for S in P_range]
end

# @load "data/compare_K_pers_rand.jld"
plt_pe = plot(framestyle=:box, grid=false, legend=:topright,
    xlabel="Max. Total Supply", ylabel="Prob. to sustain steady state", fontsize=12, size=(450, 360));

for (i, K) in pairs(K_range)
    curve = vols_all[i] ./ devols_all[i]
    idx = [findlast(==(0.0), curve); findall(curve .> 0.0)]
    xs, ys = extract_peak(P_range, curve, 0.02)
    scatter!(plt_pe, [P_range[idx[1]]], [curve[idx[1]]], mode="markers", color=colors[i], marker=:square, label="", markersize=4, markerstrokecolor=:auto)
    scatter!(plt_pe, [xs], [ys], mode="markers", color=colors[i], marker=:circle, label="", markersize=4, markerstrokecolor=:auto)
    plot!(plt_pe, P_range[idx], curve[idx], label="K=$K", lw=1.5, xaxis=:log10, color=colors[i])
end

display(plt_pe)
plot!(plt_pe, title="Another Rand σ")
# savefig(plt_pe, "figures/compare_K_pers_trop.pdf");
@save "data/compare_K_pers_rand_another.jld" K_range P_range σ vols_all devols_all plt_pe