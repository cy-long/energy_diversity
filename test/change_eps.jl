include("../src/EnerFeas.jl")
using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots, Printf
using DataFrames

function generate_problem_simple(σ::Matrix{Float64}, d::Vector{Float64}, N0::Vector{Float64})
    k = fill(0.0, size(σ, 1))
    Λ = inv(σ); P = (Λ + Λ')/2; c = -Λ * d
    return EnergyConstrProb(σ,Λ,P,c,k,d,N0,size(σ,1),0.0,:total)
end

function create_sigma(σ, n_scale, ϵ=0.25)
    # σ = rand(Normal(0,n_var), K, K)
    c = minimum(eigvals((σ + σ') / 2), init=0.0)
    sh = -c + ϵ
    σ[diagind(σ)] .+= sh
    return(n_scale * σ)
end

Random.seed!(20);
σ = rand(Normal(0,3.0), 4, 4);
d = (fill(1.0, 4)+randn(4) * 0.2) * 10.0;
N0 = fill(1.0, 4) * 5.0;
Q_range = Vector(0.01:0.1:10.0) * 10^4;
devols = [S^4 / (factorial(4) * prod(N0)) for S in Q_range];

plt = plot(framestyle=:box, grid=false, legend=:topright,
    xlabel="Q", ylabel="Pr", fontsize=12, size=(450, 360));

for ϵ in [0.1, 0.25, 0.5, 0.75, 2.0]
    σ_diss = create_sigma(σ, 20.0, ϵ);
    p = generate_problem_simple(σ_diss, d, N0);
    vols = volume_range_EFD(p, Q_range, n_sample=25000);
    plot!(plt, Q_range, vols./devols, label=@sprintf("ϵ=%.2f", ϵ), lw=1.5)
end

display(plt)


# using LinearAlgebra
# using Plots


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