"""Demonstrate the relevance of the parameters in a realistic scenario"""
"""Putting it simpler: how is the scaling of σ, d, Q, N⁰ correlated?"""

include("../src/EnerFeas.jl")
using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots
using DataFrames

function generate_problem_scale(σ::Matrix{Float64}, d::Vector{Float64}, N0::Vector{Float64})
    k = fill(0.0, size(σ, 1))
    Λ = inv(σ); P = (Λ + Λ')/2; c = -Λ * d
    return EnergyConstrProb(σ,Λ,P,c,k,d,N0,size(σ,1),0.0,:total)
end

K = 4

# 0.1, 1.0, 10.0
n_scale = 100.0

# 1.0, 10.0, 100.0
d_scale = 100.0

# 1.0, 10.0, 100.0
N0_scale = 200.0

# 1.0, 10.0, 100.0
Q_scale = 5*10.0^5

Random.seed!(1101);
σ = generate_sigma(K, false, true, n_scale, 1.0);
d = (fill(1.0, K)+randn(K) * 0.1) * d_scale;
N0 = (fill(1.0, K)+randn(K) * 0.1) * N0_scale;
p = generate_problem_scale(σ, d, N0);

Q_range = Vector(0.01:0.25:20.0) * Q_scale;
vols = volume_range_EFD(p, Q_range, n_sample=3*10^4);
devols = [S^p.K / (factorial(p.K) * prod(p.N⁰)) for S in Q_range];

plot(Q_range, vols./devols, label="1", color=:blue, lw=1.5)

Q_opt, _ = extract_peak(Q_range, vols./devols, 0.04);

q = Q_opt / (d_scale * N0_scale^2);
τ = (N0_scale * n_scale) / d_scale;

# add to a dataframe
push!(df, (n_scale, d_scale, N0_scale, Q_scale, Q_opt, q, τ))

# scatter q and τ
scatter(df.q, df.τ, xlabel="q", ylabel="τ", label="K=4", legend=:topright, title="Scaling of q and τ")

