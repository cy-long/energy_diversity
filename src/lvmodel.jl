"""
This script parameterizes energetic Lokta Volterra systems, which can be studied later
with HR sampling and computational geometry approaches.
"""

using LinearAlgebra
using Random
using Distributions

abstract type AbstractMatrixConstraintProblem end

mutable struct EnergyConstraintProblem{T} <: AbstractMatrixConstraintProblem
    σ::Matrix{T} # (energetic) interaction matrix
    Λ::Matrix{T} # inverse interaction matrix, Λ = σ⁻¹
    Q::Matrix{T} # quadratic form, Q = (Λ + Λ')/2
    c::Vector{T} # demands offset, c = -Λ * d
    m::Vector{T} # vector of bodysize
    d::Vector{T} # vector of demand
    N⁰::Vector{T} # vector of minimal biomass 
    K::Int # number of species
    S::T # supply cap, could be individual or total
    type::Symbol # specify the type of constraint
end

"""
Generate a random energy constrained ecosystem with:
K species, S supply cap, log-normal bodymass, normal interaction strength, -1/4 demands scaling, minimal biomass, 
"""
function create_energy_problem(K::Int, S::Float64, type::Symbol; seed=nothing)
    if seed !== nothing
        Random.seed!(seed)  # Set the seed for reproducibility
    end
    if type ∉ [:Individual, :Total]
        throw(ArgumentError("type must be either :Individual or :Total"))
    end
    m = rand(LogNormal(0, 0.5), K)
    σ = abs.(randn(K, K) ./m')
    d = m.^(-1/4)
    N⁰ = fill(0.1, K)
    Λ = inv(σ)
    Q = (Λ + Λ')/2
    c = -Λ * d
    return EnergyConstraintProblem(σ,Λ,Q,c,m,d,N⁰,K,S,type)
end


"""
Add minimal self-regulation to the interaction to make σ dissipative (σ+σ' positive def.)
"""
function make_problem_dissipative!(p::EnergyConstraintProblem, ϵ::Float64=1e-3)
    σ = p.σ
    m = p.m
    M_neg = Diagonal(m.^(1/2))

    # increase self-regulation while keeping network structure
    c = -min(minimum(eigvals(M_neg * ((σ + σ')/2) * M_neg)), 0)
    σ += (c+ϵ) * Diagonal(m.^(-1)) 
    Λ = inv(σ)

    p.σ, p.Λ, p.Q, p.c = σ, Λ, (Λ + Λ') / 2, -Λ * p.d
    return(minimum(eigvals(((Λ + Λ')/2)))>0)
end


"""
Calculate the baseline supply needed to sustain the ecosystem with least biomass (averaged by K)
"""
baseline_supply(p::EnergyConstraintProblem) = sum(p.d + inv(p.Λ) * p.N⁰)/p.K

"""Calculate the individual supply at state s (averaged by K)"""
individual_supply(s::Vector{Float64}, p::EnergyConstraintProblem) = dot(s, p.m)/p.K

"""Calculate the total supply at state s (weighed by N)"""
total_supply(s::Vector{Float64}, p::EnergyConstraintProblem) = transpose(s) * p.Λ * (s - p.d)