module EnerFeas

using Random, Distributions
using SpecialFunctions
using LinearAlgebra
using JuMP, Gurobi
using Polyhedra, QHull
using Plots, ProgressMeter

export EnergyConstrProb, EcosysConfig
export check_feasible_EFD, sample_EFD, volume_EFD, volume_cascade_EFD

const GRB_ENV = Gurobi.Env(output_flag=0)

# instead of "problem", these are just paramers. We need to reconfigurate this. 
# Ideally allow people to perfom whatever constraints they like to do
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

function create_poly(p::EnergyConstrProb)
    if p.type == :indiv
        A = -vcat(I, p.Λ, -p.m')
        b = -vcat(zeros(p.K), p.N⁰-p.c, -p.K*p.S)
    elseif p.type == :total # not a finite polyhedron
        A = -vcat(I, p.Λ)
        b = -vcat(zeros(p.K), p.N⁰-p.c)
    end
    return polyhedron(hrep(A, b))
end

function make_isotropic(p::EnergyConstrProb)
    if p.type == :total
        L = cholesky(p.Q).U
        invL = inv(L)
        A = -vcat(I, p.Λ) * invL
        b = -vcat(zeros(p.K), p.N⁰ - p.c)
        yc = -0.5 * invL' * p.c
        t = p.S + 0.25 * p.c' * inv(p.Q) * p.c
    elseif p.type == :indiv
        L = Matrix(1.0I, p.K, p.K) # freeze the transformation functionality; not needed and causing issues
        invL = L
        A = -vcat(I, p.Λ, -p.m') * invL
        b = -vcat(zeros(p.K), p.N⁰ - p.c, -p.K*p.S)
        yc = fill(NaN, p.K); t = NaN
    end
    return IsoTrans(L, invL, yc, t, A, b)
end

function check_feasible_EFD(p::EnergyConstrProb)
    A = -vcat(I, p.Λ)
    b = -vcat(zeros(p.K), p.N⁰ - p.c)
    model = Model()
    @variable(model, s[i = 1:p.K])

    if p.type == :indiv
        A = vcat(A, p.m')
        b = vcat(b, p.K*p.S)
    elseif p.type == :total
        @constraint(model, s'*p.Q*s + p.c'*s .<= p.S)
    end
    @constraint(model, A*s .<= b)
    @objective(model, Max, 0) # dummy objective
    set_optimizer(model, () -> Gurobi.Optimizer(GRB_ENV)); #> handling usage outside this package
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "LogToConsole", 0)
    optimize!(model)
    if termination_status(model) != MOI.OPTIMAL
        return false
    end

    @objective(model, Max, sum(s))
    optimize!(model)
    if termination_status(model) == MOI.OPTIMAL
        return objective_value(model) > 1e-9
    else
        return false
    end
end


# go_back only for testing purposes
function sample_EFD(p::EnergyConstrProb, N::Int=10000, go_back::Bool=true)
    check_feasible_EFD(p) || throw(ErrorException("Cannot sample on infeasible domain"))
    itp = make_isotropic(p)
    if p.type == :indiv
        domain = InterPolySpheres(itp.A, itp.b, [], :indiv)
    elseif p.type == :total
        sp0 = Sphere(itp.yc, sqrt(itp.t))
        domain = InterPolySpheres(itp.A, itp.b, [sp0], :total)
    end
    nt = 10; ns = ceil(Int, 2*N/nt) #TODO: find the best nt and ns
    samples = hr_sample(domain, nt, ns, chevball(domain))
    if go_back 
        samples = [itp.invL * y for y in samples] # transform back to s space
    end
    if (length(samples) > N)
        samples = samples[1:N] # truncate to N samples
    end
    return samples
end


# higher level wrapper for EnergyConstrProb, just specify the domains
function volume_EFD(p::EnergyConstrProb, exact::Bool=false, N::Int=10, eN::Int=1)
    check_feasible_EFD(p) || return 0.0, NaN

    itp = make_isotropic(p)

    if p.type == :indiv
        domain = InterPolySpheres(itp.A, itp.b, [], :indiv)
        vol_full = (p.K*p.S)^p.K / (factorial(p.K) * prod(p.m))

    elseif p.type == :total
        sp0 = Sphere(itp.yc, sqrt(itp.t))
        domain = InterPolySpheres(itp.A, itp.b, [sp0], :total)
        exact && throw(ErrorException("Exact volume not supported for Total Energy Constraint"))
        vol_full = vol_sphere(sp0) * det(itp.invL)
    end

    vol_in_y = volume_domain(domain, N, eN, exact)
    vol_in_s = vol_in_y * det(itp.invL)
    return vol_in_s, vol_full
end


# p0: domain of the last (biggest) region, S_range: increasing range of total supply values
function volume_cascade_EFD(p0::EnergyConstrProb, S_range::Vector{Float64}, α::Float64=0.5)
    p0.type == :indiv && throw(ErrorException("Cascade volume not supported for Individual Energy Constraint"))
    S_range[1] > S_range[end] || reverse!(S_range)
    
    itp = make_isotropic(p0)
    quadbounds = [Sphere(itp.yc, sqrt(itp.t + S - p0.S)) for S in S_range]
    regions = [InterPolySpheres(itp.A, itp.b, [q], :total) for q in quadbounds]
    balls = [chevball(r) for r in regions]
    
    volumes = [0.0 for _ in S_range]
    volumes[1] = volume_domain(regions[1], 10, 1, true); r = 1.0
    
    @showprogress desc="Volume: $(S_range[end]) to $(S_range[1])" for i in eachindex(regions)
        if r < 1e-6 || i == length(regions) || balls[i].r < 1e-9
            break
        end
        samples_i = hr_sample(regions[i], 10, 5000, balls[i])
        inside_R_ii = [is_inside(s, regions[i+1]) for s in samples_i]
        r = mean(inside_R_ii)
        b = mean([is_inside_sphere(s, balls[i+1]) for s in samples_i])
        volumes[i+1] = α * (volumes[i] * r) + (1 - α) * (vol_sphere(balls[i+1]) * r / b)
    end
    println("")

    reverse!(S_range); reverse!(volumes)
    volumes[isnan.(volumes)] .= 0.0
    return volumes * det(itp.invL) # transform back to s space
end

include("core.jl")
include("generate.jl")
include("visualize.jl")


# temporarily for development!!
export IsoTrans, create_poly, make_isotropic, Sphere, rand_sphere, vol_sphere, InterPolySpheres, chevball, hr_step, hr_sample, is_inside, is_inside_sphere, volume_domain, show_chevball, show_linear, show_quadratic, grow_quadratic

export ecosys_config, generate_sigma_arrays, generate_problem, baseline_supply

end