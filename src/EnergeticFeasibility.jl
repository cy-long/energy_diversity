module EnergeticFeasibility

using Random, Distributions
using SpecialFunctions
using LinearAlgebra
using JuMP, Gurobi
using Polyhedra, QHull, MinimumVolumeEllipsoids
using Plots

export EnergyConstrProb
export sample_EFD, volume_EFD
export show_chevball, show_linear

const GRB_ENV = Gurobi.Env(output_flag=0)

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

# linear transformed model using L; variables are defined at y instead of s space also.
struct IsoTrans
    L::Matrix{Float64} # Cholesky upper matrix, Q=L'L
    invL::Matrix{Float64} # inverse of L
    yc::Vector{Float64} # transformed center of sphere
    t::Float64 # square of radius of sphere
    A::Matrix{Float64} # transformed linear constraints
    b::Vector{Float64} # linear constraints constants
end

#actually the interface to translate p into geometrical stuff
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

# go_back only for testing purposes
function sample_EFD(p::EnergyConstrProb, N::Int=10000, go_back::Bool=true)
    itp = make_isotropic(p)
    if p.type == :Individual
        domain = InterPolySpheres(itp.A, itp.b, [], :Individual)
    elseif p.type == :Total
        sp0 = Sphere(itp.yc, sqrt(itp.t))
        domain = InterPolySpheres(itp.A, itp.b, [sp0], :Total)
    end
    nt = 10; ns = ceil(Int, 2*N/nt) #TODO: find the best nt and ns
    samples = hr_sample(domain, nt, ns, chevball(domain))
    if go_back 
        samples = [itp.invL * y for y in samples] # transform back to s space
    end
    if (length(samples) > N)
        samples = samples[1:N] # truncate to N samples
    end
    return(samples)
end


# higher level wrapper for EnergyConstrProb, just specify the domains
function volume_EFD(p::EnergyConstrProb, N::Int=10, eN::Int=1, exact::Bool=false)
    itp = make_isotropic(p)

    if p.type == :Individual
        domain = InterPolySpheres(itp.A, itp.b, [], :Individual)
    elseif p.type == :Total
        sp0 = Sphere(itp.yc, sqrt(itp.t))
        domain = InterPolySpheres(itp.A, itp.b, [sp0], :Total)
        exact && throw(ErrorException("Exact volume not supported for Total Energy Constraint"))
    end

    vol_in_y = volume_domain(domain, N, eN, exact)
    vol_in_s = vol_in_y * det(itp.invL)
    return(vol_in_s)
end

# add other functions that study the property of EnergyConstrProb, for instance, when to have feasible solution

include("core.jl")
include("visualization.jl")

end