""" testing """

include("src/EnergeticFeasibility.jl")
using .EnergeticFeasibility
using Random, Distributions, LinearAlgebra, Plots


function test_quad(K::Int=3, S::Float64=4.0, ϵ::Float64=0.75, seed::Int=42)
    Random.seed!(seed)
    σ = randn(K,K)
    # σ = Matrix(1.0I, K, K)
    c = minimum(eigvals((σ + σ')/2))
    if c < 0
        σ += (-c + ϵ) * Diagonal(fill(1,K))
    end
    Λ = inv(σ)
    Q = (Λ + Λ')/2
    m = fill(1.0,K); d = fill(1.0, K); N⁰ = fill(1e-3, K); c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:Total)
end

function test_linear(K::Int=3, S::Float64=4.0, ϵ::Float64=0.75, seed::Int=42)
    Random.seed!(seed)
    σ = randn(K,K)
    c = minimum(eigvals((σ + σ')/2))
    if c < 0
        σ += (-c + ϵ) * Diagonal(fill(1,K))
    end
    Λ = inv(σ)
    Q = (Λ + Λ')/2
    m = fill(1.0,K); d = fill(1.0, K); N⁰ = fill(1e-3, K); c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:Individual)
end


#1 sampling linear constraint domain
p = test_linear(2, 2.5, 0.5, 32);
samples = sample_EFD(p, 15000);

plot(ratio=1)
show_linear(p, samples)
show_chevball(p)

#2 sampling quadratic constraint domain
p = test_quad(2, 20.0, 0.75, 42);
samples = sample_EFD(p, 15000);

plot(ratio=1)
show_quadratic(p, samples)

#3 linear volume estimation
p = test_linear(5, 7.5, 0.5, 4);
volume_EFD(p, 10, 1, true)
volume_EFD(p, 10, 1, false)

#4 quadratic volume estimation
p1 = test_quad(4, 8.5, 0.75, 42);
volume_EFD(p1, 10, 1, true)
volume_EFD(p1, 10, 1, false)