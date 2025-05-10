module EnerFeas

using Random, Distributions
using SpecialFunctions
using LinearAlgebra
using JuMP, Gurobi
using Polyhedra, QHull
using Plots, ProgressMeter

const GRB_ENV = Gurobi.Env(output_flag=0)

# instead of "problem", these are just paramers. We need to reconfigurate this. 
# Ideally allow people to perfom whatever constraints they like to do
mutable struct EnergyConstrProb{T}
    σ::Matrix{T} # (energetic) interaction matrix
    Λ::Matrix{T} # inverse interaction matrix, Λ = σ⁻¹
    Q::Matrix{T} # quadratic form, Q = (Λ + Λ')/2
    c::Vector{T} # demands offset, c = -Λ * d
    k::Vector{T} # vector of coefficients for feasibility tolerance
    d::Vector{T} # vector of demand
    N⁰::Vector{T} # vector of minimal biomass 
    K::Int # number of species
    S::T # supply cap, could be individual or total
    type::Symbol # specify the type of constraint
end

# linear transformed model using L; variables are defined at y instead of s space also.
struct IsotropicTransParams
    L::Matrix{Float64} # Cholesky upper matrix, Q=L'L
    invL::Matrix{Float64} # inverse of L
    yc::Vector{Float64} # transformed center of sphere
    t::Float64 # square of radius of sphere
    A::Matrix{Float64} # transformed linear constraints
    b::Vector{Float64} # linear constraints constants
end

function translate_EFD(p::EnergyConstrProb)
    L = Matrix(1.0I, (p.K, p.K)); invL = Matrix(1.0I, (p.K, p.K))
    yc = fill(NaN, p.K); t = NaN
    if p.type == :total
        L = cholesky(p.Q).U
        invL = inv(L)
        yc = -0.5 * invL' * p.c
        t = p.S + 0.25 * p.c' * inv(p.Q) * p.c
        A = -vcat(I, p.Λ, -p.N⁰') * invL
        b = -vcat(zeros(p.K), p.k .* p.N⁰ - p.c, -p.S)
    elseif p.type == :indiv
        A = -vcat(I, p.Λ, -p.N⁰')
        b = -vcat(zeros(p.K), p.k .* p.N⁰ - p.c, -p.S)
    end
    return IsotropicTransParams(L, invL, yc, t, A, b)
end

function create_poly(p::EnergyConstrProb)
    if p.type == :indiv
        A = -vcat(I, p.Λ, -p.N⁰)
        b = -vcat(zeros(p.K), p.k .* p.N⁰-p.c, -p.S)
    elseif p.type == :total # not a finite polyhedron
        A = -vcat(I, p.Λ)
        b = -vcat(zeros(p.K), p.k .* p.N⁰-p.c)
    end
    return polyhedron(hrep(A, b))
end

function check_feasible_EFD(p::EnergyConstrProb)
    model = Model()
    @variable(model, s[i = 1:p.K])

    A = -vcat(I, p.Λ, -p.N⁰')
    b = -vcat(zeros(p.K), p.k .* p.N⁰ - p.c, -p.S)
    @constraint(model, A*s .<= b)

    if p.type == :total
        @constraint(model, s'*p.Q*s + p.c'*s <= p.S)
    end

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
    itp = translate_EFD(p)
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
    check_feasible_EFD(p) || return 0.0

    itp = translate_EFD(p)

    if p.type == :indiv
        domain = InterPolySpheres(itp.A, itp.b, [], :indiv)

    elseif p.type == :total
        sp0 = Sphere(itp.yc, sqrt(itp.t))
        domain = InterPolySpheres(itp.A, itp.b, [sp0], :total)
        exact && throw(ErrorException("Exact volume not supported for Total Energy Constraint"))
        # vol_full = vol_sphere(sp0) * det(itp.invL)
    end

    vol_in_y = volume_domain(domain, N, eN, exact)
    vol_in_s = vol_in_y * det(itp.invL)
    # return vol_in_s, vol_full
    return vol_in_s
end

# p: a seed problem specifying the unchanged constraints; P_range: range of total supply that changes
# n_thread: inside each hr_sample, # thread of each sampling [future: change to parallel]
# n_sample: inside each hr_sample, # samples for inside each thread
# α: [experimental] balance parameter from cascade and ball vol estimation
function volume_range_EFD(p::EnergyConstrProb, P_range::Vector{Float64}; n_thread::Int=10, n_sample::Int=10^4, n_layer::Int=10, α::Float64=1.0)
    if p.type == :indiv 
        volumes = [0.0 for _ in P_range]
        @showprogress desc="Volume: $(P_range[1]) to $(P_range[end])" for (i,S) in pairs(P_range)
            p.S = S
            volumes[i] = volume_EFD(p, true)
        end

    elseif p.type == :total
        volumes = [0.0 for _ in P_range]
        P_range[1] > P_range[end] || reverse!(P_range)
        
        itp = translate_EFD(p)
        quadbounds = [Sphere(itp.yc, sqrt(S + itp.t-p.S))
            for S in P_range] 
        regions = [InterPolySpheres(itp.A, vcat(itp.b[1:end-1], itp.b[end]+S-p.S), [q], :total)
            for (S, q) in zip(P_range, quadbounds)]
        balls = [chevball(r) for r in regions]

        volumes[1] = volume_domain(regions[1], n_layer, 1, true)
        if volumes[1] == 0.0
            return volumes
        end

        r = 1.0;
        @showprogress desc="Volume: $(P_range[end]) to $(P_range[1])" for i in eachindex(regions)
            if r < 1e-6 || i == length(regions) || balls[i].r < 1e-9
                break
            end
            samples_i = hr_sample(regions[i], n_thread, n_sample, balls[i])
            inside_R_ii = [is_inside(s, regions[i+1]) for s in samples_i]
            r = mean(inside_R_ii)
            b = mean([is_inside_sphere(s, balls[i+1]) for s in samples_i])
            volumes[i+1] = α * (volumes[i] * r) + (1 - α) * (vol_sphere(balls[i+1]) * r / b)
        end

        reverse!(P_range); reverse!(volumes)

        volumes[isnan.(volumes)] .= 0.0
        volumes = volumes * det(itp.invL) # back in s-space
    end

    return volumes
end

# Calculate the baseline supply needed to sustain the ecosystem with least biomass
baseline_supply(p::EnergyConstrProb) = dot(p.d + p.σ * (p.k .* p.N⁰), p.N⁰)

# Calculate the total supply at state s (weighed by N)
total_supply(s::Vector{Float64}, p::EnergyConstrProb) = transpose(s) * p.Λ * (s - p.d)

include("core.jl")
include("generate.jl")
include("visualize.jl")
include("analyze.jl")

export EnergyConstrProb, EcosysConfig
export check_feasible_EFD, sample_EFD, volume_EFD, volume_range_EFD
export baseline_supply
export ecosys_config, generate_sigma_arrays, generate_problem, generate_sigma, generate_trophic
export extract_critical, extract_peak, score_saturation, score_unimodal

# ---- these exports for temporary development only ----
export show_linear, show_quadratic
export IsotropicTransParams, create_poly, translate_EFD
export Sphere, rand_sphere, vol_sphere, InterPolySpheres, chevball, hr_step, hr_sample, is_inside, is_inside_sphere, volume_domain, show_chevball, grow_quadratic
export smooth, smooth_curve
# ---- ----

end