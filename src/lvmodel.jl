using LinearAlgebra
using Random
using Distributions

abstract type AbstractMatrixConstraintProblem end

struct EnergyConstraintProblem{T} <: AbstractMatrixConstraintProblem
    Λ::Matrix{T} # inverse interaction matrix
    N⁰::Vector{T} # baseline biomass vector
    m::Vector{T} # bodysize vector
    K::Int # number of species
    d::Vector{T} # demand vector
    Sᵢ::T # average individual supply cap
end

mutable struct HRSampler{T <: AbstractMatrixConstraintProblem}
    # Problem information
    pblm::T
    Z::Matrix{Float64} # nullspace of pblm.E
    warmup::Vector{Vector{Float64}}
    center::Vector{Float64}
    n_warmup::Int
    n_samples::Int
    prev::Vector{Float64}
    feas_tol::Float64
    bounds_tol::Float64
    feasible::Bool
    RNG::Random.MersenneTwister
end

function generate_inte_matrix(Ks::Int)
    inte_matrix = -abs.(randn(Ks, Ks))
    inte_matrix[diagind(inte_matrix)] .= -1
    return inte_matrix
end

function create_energy_problem(K::Int, Sind::Float64; seed=nothing)
    if seed !== nothing
        Random.seed!(seed)  # Set the seed for reproducibility
    end
    m = rand(LogNormal(0, 0.5), K)
    σ = -generate_inte_matrix(K) ./m' # currently purely compatitive
    Λ = inv(σ)
    d = m.^(-1/4)
    N⁰ = fill(1e-2, K)
    Sᵢ = Sind
    # m = fill(1.0, K)
    # N⁰ = fill(1e-2, K)
    # d = fill(0.0, K)
    return EnergyConstraintProblem(Λ, N⁰, m, K, d, Sᵢ)
end

function create_sampler(problem::EnergyConstraintProblem)
    pblm = problem
    Z = Matrix(Diagonal(ones(pblm.K))) # no equality const., null space=the entire space
    warmup = Vector{Vector{Float64}}()
    center = zeros(pblm.K)
    n_warmup = 0
    n_samples = 0
    prev = zeros(pblm.K)
    feas_tol = 1e-6
    bounds_tol = 1e-6
    feasible = false
    RNG = MersenneTwister(42)

    return HRSampler(pblm, Z, warmup, center, n_warmup, n_samples, prev, feas_tol, bounds_tol, feasible, RNG)
end

function baseline_supply(problem::EnergyConstraintProblem)
    return sum(problem.d + inv(problem.Λ) * problem.N⁰)/problem.K
end