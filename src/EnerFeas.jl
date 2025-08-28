module EnerFeas

using Random, Distributions
using SpecialFunctions
using LinearAlgebra
using JuMP, SCS
using Polyhedra, QHull
using Plots
using Printf, ProgressMeter
using Loess

# instead of "problem", these are just paramers. We need to reconfigurate this. 
# Ideally allow people to perfom whatever constraints they like to do
mutable struct EnergyConstrProb{T}
    σ::Matrix{T} # (energetic) interaction matrix
    Λ::Matrix{T} # inverse interaction matrix, Λ = σ⁻¹
    Q::Matrix{T} # quadratic form, Q = (Λ + Λ')/2
    c::Vector{T} # demands offset, c = -Λ * d
    ϵ::Vector{T} # coefficients of threshold in steady state biomass
    d::Vector{T} # vector of demand
    N⁰::Vector{T} # vector of minimal viable biomass 
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
        b = -vcat(zeros(p.K), p.ϵ .* p.N⁰ - p.c, -p.S)
    elseif p.type == :indiv
        A = -vcat(I, p.Λ, -p.N⁰')
        b = -vcat(zeros(p.K), p.ϵ .* p.N⁰ - p.c, -p.S)
    end
    return IsotropicTransParams(L, invL, yc, t, A, b)
end

function create_poly(p::EnergyConstrProb)
    if p.type == :indiv
        A = -vcat(I, p.Λ, -p.N⁰)
        b = -vcat(zeros(p.K), p.ϵ .* p.N⁰-p.c, -p.S)
    elseif p.type == :total # not a finite polyhedron
        A = -vcat(I, p.Λ)
        b = -vcat(zeros(p.K), p.ϵ .* p.N⁰-p.c)
    end
    return polyhedron(hrep(A, b))
end

function check_feasible_EFD(p::EnergyConstrProb)
    model = Model()
    @variable(model, s[i = 1:p.K])

    A = -vcat(I, p.Λ, -p.N⁰')
    b = -vcat(zeros(p.K), p.ϵ .* p.N⁰ - p.c, -p.S)
    @constraint(model, A*s .<= b)

    if p.type == :total
        @constraint(model, s'*p.Q*s + p.c'*s <= p.S)
    end

    @objective(model, Max, 0) # dummy objective
    set_optimizer(model, SCS.Optimizer);
    set_optimizer_attribute(model, "verbose", 0)
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
        domain = InterPolyBalls(itp.A, itp.b, [], :indiv)
    elseif p.type == :total
        sp0 = Sphere(itp.yc, sqrt(itp.t))
        domain = InterPolyBalls(itp.A, itp.b, [sp0], :total)
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
        domain = InterPolyBalls(itp.A, itp.b, [], :indiv)

    elseif p.type == :total
        sp0 = Sphere(itp.yc, sqrt(itp.t))
        domain = InterPolyBalls(itp.A, itp.b, [sp0], :total)
        exact && throw(ErrorException("Exact volume not supported for Total Energy Constraint"))
        # vol_full = vol_sphere(sp0) * det(itp.invL)
    end

    vol_in_y = volume_domain(domain, N, eN, exact)
    vol_in_s = vol_in_y * det(itp.invL)
    # return vol_in_s, vol_full
    return vol_in_s
end

# p: a seed problem specifying the unchanged constraints; Q_range: range of total supply that changes
# n_thread: number thread of each sampling [future: make it parallel]
# n_sample: number of sample per thread
# n_layer: number of estimation phases per volume calculation
function volume_range_EFD(p::EnergyConstrProb, Q_range::Vector{Float64}; n_thread::Int=10, n_sample::Int=10^4, n_layer::Int=10, show_p::Bool=true, show_dt::Float64=0.5)
    desc = @sprintf("Volume: %.2f to %.2f", Q_range[1], Q_range[end])
    prog = show_p ? Progress(length(Q_range), dt = show_dt, desc=desc, showspeed=true) : nothing
    
    if p.type == :indiv 
        volumes = [0.0 for _ in Q_range]
        for (i,S) in pairs(Q_range)
            p.S = S
            volumes[i] = volume_EFD(p, true)
            if prog !== nothing
                next!(prog)
            end
        end

    elseif p.type == :total
        volumes = [0.0 for _ in Q_range]
        Q_range[1] > Q_range[end] || reverse!(Q_range)
        
        itp = translate_EFD(p)
        quadbounds = [Sphere(itp.yc, sqrt(S + itp.t-p.S))
            for S in Q_range] 
        regions = [InterPolyBalls(itp.A, vcat(itp.b[1:end-1], itp.b[end]+S-p.S), [q], :total)
            for (S, q) in zip(Q_range, quadbounds)]
        balls = [chevball(r) for r in regions]

        volumes[1] = volume_domain(regions[1], n_layer, 1, true)
        if volumes[1] == 0.0
            reverse!(Q_range)
            return volumes
        end

        r = 1.0;
        for i in eachindex(regions)
            if i == length(regions) || balls[i].r < 5*1e-5 || r < 1e-6
                break
            end
            
            samples_i = nothing
            try
                samples_i = hr_sample(regions[i], n_thread, n_sample, balls[i])
            catch e
                @info "Sampling stalled at region $i, cutoff as zero volume."
                break
            end

            r = mean([is_inside(s, regions[i+1]) for s in samples_i])
            volumes[i+1] = volumes[i] * r
            if prog !== nothing
                next!(prog)
            end
        end
        if prog !== nothing
            finish!(prog)
        end

        reverse!(Q_range); reverse!(volumes)

        volumes[isnan.(volumes)] .= 0.0
        volumes = volumes * det(itp.invL) # back in s-space
    end
    return volumes
end

include("core.jl")
include("generate.jl")
include("analyze.jl")

export EnergyConstrProb, EcosysConfig
export check_feasible_EFD, sample_EFD, volume_EFD, volume_range_EFD, translate_EFD
export baseline_supply, optimal_supply
export ecosys_config, generate_sigma_arrays, generate_problem
export extract_critical, extract_peak, score_saturation, score_unimodal
export generate_model_system, sub_model_system, select_range, volume_range_flux

# ---- these lower-end exports for temporary development only ----
export Sphere, IsotropicTransParams, InterPolyBalls
export chevball, volume_domain, vol_sphere
export smooth, smooth_curve, generate_sigma
export critical_energy
export rand_sphere, hr_step, hr_sample, is_inside
export rand_warmup
# ---- ----

end