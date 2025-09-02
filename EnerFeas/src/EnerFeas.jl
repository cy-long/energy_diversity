module EnerFeas

using Random, Distributions
using SpecialFunctions
using LinearAlgebra
using JuMP, SCS
using Polyhedra, QHull
using Plots
using Printf, ProgressMeter
using Loess

mutable struct EcosysParams # (p)
    S::Int # number of populations
    σ::Matrix{Float64} # energy exchange rate matrix
    d::Vector{Float64} # energy demand rate vector
    Q::Float64 # total energy supply
    N⁰::Vector{Float64} # minimal biomass vector
    ϵ::Vector{Float64} # threshold biomass ratio vector
end

struct TransParams # (tp)
    S::Int # number of populations
    Q::Float64 # total energy supply
    Λ::Matrix{Float64} # inverse of σ
    P::Matrix{Float64} # symmetrization of Λ, P=(Λ+Λ')/2
    c::Vector{Float64} # demands offset, c = -Λ * d
    cholP::Union{Cholesky{Float64, Matrix{Float64}}, Nothing} # Cholesky factorization of P or nothing
    detL::Float64 # determinant of L, where P=L'L
    yc::Vector{Float64} # transformed center of sphere
    t::Float64 # square of radius of sphere
    A::Matrix{Float64} # transformed linear constraints
    b::Vector{Float64} # linear constraints constants
    type::Symbol # specify the type of constraint (:init, :matr)
end

function translate_EFD(p::EcosysParams, type::Symbol)::TransParams
    Λ = inv(p.σ)
    P = 0.5 * (Λ + Λ')
    c = -Λ * p.d
    if type == :matr
        if minimum(eigvals(P)) < 1e-12
            throw(ArgumentError("The system is not dissipative, check σ"))
        end
        cholP = cholesky(P)
        yc = -0.5 * (c' / cholP.U)'
        t = p.Q + 0.25 * c' * (P \ c)
        detL = prod(diag(cholP.U))
        A = -vcat(I, Λ, -p.N⁰') / cholP.U
        b = -vcat(zeros(p.S), p.ϵ .* p.N⁰ - c, -p.Q)
    else
        cholP = nothing
        yc = fill(NaN, p.S)
        t = NaN
        detL = 1.0
        A = -vcat(I, Λ, -p.N⁰')
        b = -vcat(zeros(p.S), p.ϵ .* p.N⁰ - c, -p.Q)    
    end

    return TransParams(p.S, p.Q, Λ, P, c, cholP, detL, yc, t, A, b, type)
end

function translate_domain(tp::TransParams)::InterPolyBalls
    if tp.type == :init
        return InterPolyBalls(tp.A, tp.b, [], :linr)
    elseif tp.type == :matr
        sp0 = Sphere(tp.yc, sqrt(tp.t))
        return InterPolyBalls(tp.A, tp.b, [sp0], :quad)
    else
        throw(ArgumentError("Unknown TransParams type: $(tp.type)"))
    end
end

function check_feasible_EFD(tp::TransParams)::Bool
    return is_feasible(translate_domain(tp))
end

function sample_EFD(tp::TransParams; n_sample::Int=10000, n_thread::Int=1)
    check_feasible_EFD(tp) || throw(ErrorException("The ecosystem has infeasible domain"))
    domain = translate_domain(tp)

    n_spt = ceil(Int, 2 * n_sample / n_thread)
    samples = hr_sample(domain, n_thread, n_spt, chevball(domain))

    if (length(samples) > n_sample)
        samples = samples[1:n_sample] # truncate to n_sample
    end
    if tp.type == :matr
        samples = [tp.cholP.U \ s for s in samples]# back to s-space
    end
    return samples
end


function volume_EFD(tp::TransParams; exact::Bool=false, N::Int=10, eN::Int=1)
    check_feasible_EFD(tp) || return 0.0
    if tp.type == :matr && exact
        throw(ErrorException("Exact volume not supported for Total Energy Constraint"))
    end
    domain = translate_domain(tp)
    vol_in_y = volume_domain(domain, N, eN, exact)
    vol_in_s = vol_in_y / tp.detL
    return vol_in_s
end

function volume_range_EFD(
    p::EcosysParams, type::Symbol, Q_range::Vector{Float64};
    n_thread::Int=10, n_sample::Int=10^4, n_layer::Int=10, show_prog::Bool=true, prog_dt::Float64=0.5)

    desc = @sprintf("Volume: %.2f to %.2f", Q_range[1], Q_range[end])
    prog = show_prog ? Progress(length(Q_range), dt = prog_dt, desc=desc, showspeed=true) : nothing
    
    if type == :init 
        volumes = [0.0 for _ in Q_range]
        
        for (i,Q) in pairs(Q_range)
            p.Q = Q
            volumes[i] = volume_EFD(p, :init, exact=true)
            if prog !== nothing
                next!(prog)
            end
        end

    elseif type == :matr
        tp = translate_EFD(p, :matr)
        t0 = tp.t - p.Q
        volumes = [0.0 for _ in Q_range]
        Q_range[1] > Q_range[end] || reverse!(Q_range)

        quadbounds = [Sphere(tp.yc, sqrt(Q + t0)) for Q in Q_range] 
        regions = [InterPolyBalls(tp.A, vcat(tp.b[1:end-1], tp.b[end]+Q-p.Q), [q], :quad)
            for (Q, q) in zip(Q_range, quadbounds)]
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
        volumes = volumes / tp.detL # back in s-space
    end
    return volumes
end

function volume_range_C(p::EcosysParams, Q_range::Vector{Float64})
    return [Q^p.S/(prod(p.N⁰)*factorial(p.S)) for Q in Q_range]
end

check_feasible_EFD(p::EcosysParams, type::Symbol) = check_feasible_EFD(translate_EFD(p, type))
sample_EFD(p::EcosysParams, type::Symbol; kwargs...) = sample_EFD(translate_EFD(p, type); kwargs...)
volume_EFD(p::EcosysParams, type::Symbol; kwargs...) = volume_EFD(translate_EFD(p, type); kwargs...)

include("core.jl")
include("generate.jl")
include("analyze.jl")

export EcosysParams, TransParams
export translate_EFD, translate_domain
export check_feasible_EFD, sample_EFD, volume_EFD, volume_range_EFD, volume_range_C
export generate_model_system, sub_model_system, select_range, ecosys_config, generate_sigma_arrays
export extract_critical, extract_peak, score_saturation, score_unimodal
# export Sphere, InterPolyBalls, chevball, volume_domain
export optimal_supply, baseline_supply
export InterPolyBalls, Sphere, volume_domain

end