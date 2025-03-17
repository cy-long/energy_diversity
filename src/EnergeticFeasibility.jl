module EnergyFeasibility

using Random, Distributions
using SpecialFunctions
using LinearAlgebra
using JuMP, Gurobi
using Polyhedra, QHull, MinimumVolumeEllipsoids
using Plots

export EnergyConstrProb, IsoTrans
export create_poly, volume_EFD


mutable struct EnergyConstrProb{T}
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


# linear transformed model using L; variables are defined at y instead of s space also.
struct IsoTrans
    L::Matrix{Float64} # Cholesky upper matrix, Q=L'L
    invL::Matrix{Float64} # inverse of L
    yc::Vector{Float64} # transformed center of sphere
    t::Float64 # square of radius of sphere
    A::Matrix{Float64} # transformed linear constraints
    b::Vector{Float64} # linear constraints constants
end

struct Sphere
    c::Vector{Float64} # center
    r::Float64 # radius
end

mutable struct InterPolySpheres
    A::Matrix{Float64}
    b::Vector{Float64}
    sps::Vector{Sphere}
    chev::Sphere
    type::Symbol # :Individual or :Total
end

function vol_sphere(sp::Sphere)
    K = length(sp.c)
    return (pi^(K/2) / gamma(K/2 + 1)) * sp.r^K
end

function chevball(p::EnergyConstrProb, itp::IsoTrans)
    A = itp.A; b = itp.b
    A_norms = [norm(row) for row in eachrow(A)]

    model = Model()
    @variable(model, x[i = 1:p.K])
    @variable(model, r >= 0)
    @constraint(model, A * x + A_norms * r .<= b)
    if p.type == :Total
        yc = itp.yc;  t = itp.t
        @constraint(model, [sqrt(t)-r; x - yc] in SecondOrderCone())
    end
    @objective(model, Max, r)

    set_optimizer(model, () -> Gurobi.Optimizer());
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "LogToConsole", 0)
    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        error("Chevball optimization failed: $(termination_status(model))")
    else
        return Sphere(value.(x), value.(r))
    end
end

function min_vol_ellipsoid(po::Polyhedron)
    pts = hcat(points(po)...)
    ϵ = minimum_volume_ellipsoid(pts)
    return Matrix(ϵ.H)
end

function create_poly(p::EnergyConstrProb)
    if p.type == :Individual
        A = -vcat(I, p.Λ, -p.m')
        b = -vcat(zeros(p.K), p.N⁰-p.c, -p.K*p.S)
    elseif p.type == :Total # not a finite polyhedron
        A = -vcat(I, p.Λ)
        b = -vcat(zeros(p.K), p.N⁰-p.c)
    end
    return(polyhedron(hrep(A, b)))
end

function make_isotropic(p::EnergyConstrProb)
    if p.type == :Total
        L = cholesky(p.Q).U
        invL = inv(L)
        A = -vcat(I, p.Λ) * invL
        b = -vcat(zeros(p.K), p.N⁰ - p.c)
        yc = -0.5 * invL' * p.c
        t = p.S + 0.25 * p.c' * inv(p.Q) * p.c
    elseif p.type == :Individual
        po = create_poly(p)
        H = min_vol_ellipsoid(po) # do we have to do this?
        L = cholesky(H).U
        invL = inv(L)
        A = -vcat(I, p.Λ, -p.m') * invL
        b = -vcat(zeros(p.K), p.N⁰ - p.c, -p.K*p.S)
        yc = fill(NaN, p.K); t = NaN
    end
    return IsoTrans(L, invL, yc, t, A, b)
end

include("sampler.jl")

end