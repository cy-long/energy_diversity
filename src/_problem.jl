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