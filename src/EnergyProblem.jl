module EnergyProblem

export EnergyConstraintProblem, HRSampler, IsotropicTransformParams

using Random
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

mutable struct HRSampler{T <: AbstractMatrixConstraintProblem}
    # Problem information
    problem::T
    warmup::Vector{Vector{Float64}}
    start::Vector{Float64}
    prev::Vector{Float64}
    n_warmup::Int
    n_samples::Int
    abs_tol::Float64
    scale_tol::Float64
    feasible::Bool
    RNG::Random.MersenneTwister
end

# linear transformed model using L; variables are defined at y instead of s space also.
struct IsotropicTransformParams
    L::Matrix{Float64} # Cholesky upper matrix, Q=L'L
    invL::Matrix{Float64} # inverse of L
    yc::Vector{Float64} # transformed center of sphere
    t::Float64 # square of radius of sphere
    A::Matrix{Float64} # transformed linear constraints
    b::Vector{Float64} # linear constraints constants
end

end