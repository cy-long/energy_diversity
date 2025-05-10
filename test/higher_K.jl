include("../src/EnerFeas.jl");

using .EnerFeas
using Random, Distributions, LinearAlgebra
using ProgressMeter, Plots, SpecialFunctions, JLD2

## --- demonstrate a single case with K=8 ---

K = 8; P_range = Vector(0.01:2.0:250.0);
ec = ecosys_config(K=K, S_type=:total, k_param=0.01, d_param=1.0, N0_param=3.0, seed=100);
σ = generate_sigma_arrays(ec, 1);
p = generate_problem(ec, σ);
devols = [S^p.K / (factorial(p.K) * prod(p.N⁰)) for S in P_range];

Random.seed!(1001);
# vols0 = volume_range_EFD(p, P_range, n_sample=10^5, α=0.5) #lead to large fluctuations
# vols1 = volume_range_EFD(p, P_range, n_sample=10^5) #already quite smooth and accurate

@load "higher_K.jld" K P_range ec σ p vols0 vols1 devols
plot(ylim=(0,0.25),framestyle=:box,grid=false)
plot!(P_range, vols0./devols, label="reproduce noise", color=:grey,lw=1.5)
plot!(P_range, vols1./devols, label="accurate sampling", color=:blue,lw=1.5)

score_saturation(P_range, vols1./devols)
score_unimodal(P_range, vols1./devols)
display(current())