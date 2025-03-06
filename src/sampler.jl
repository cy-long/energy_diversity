using LinearAlgebra
using Random
using Distributions
using RecursiveArrayTools
using JuMP
using Gurobi

# replace with your Gurobi license path
ENV["GRB_LICENSE_FILE"] = "/Users/longcy/Documents/Gurobi/gurobi.lic" 

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


function create_sampler(problem::EnergyConstraintProblem; samp_seed::Int=42, abs_tol::Float64=1e-6, scale_tol::Float64=1e-6)
    warmup = Vector{Vector{Float64}}()
    start = zeros(problem.K)
    prev = zeros(problem.K)
    n_warmup = 0
    n_samples = 0
    feasible = false
    RNG = MersenneTwister(samp_seed)
    return HRSampler(problem, warmup, start, prev, n_warmup, n_samples, abs_tol, scale_tol, feasible, RNG)
end


function warmup!(samp::HRSampler)
    p = samp.problem
    abs_tol= samp.abs_tol
    scale_tol = samp.scale_tol
    N_eff = p.N⁰ - p.c
    
    # Initialize a JuMP Model with Gurobi solver for the warmup
    model = Model(() -> Gurobi.Optimizer(GRB_ENV));
    set_optimizer_attribute(model, "OutputFlag", 0)   # Mute all Gurobi output
    set_optimizer_attribute(model, "LogToConsole", 0) # Prevent logging to console

    # Specify the variables and constraints for the warmup
    @variable(model, s[i = 1:p.K] >= abs_tol);
    @constraint(model, p.Λ * s .>= N_eff * (1 + scale_tol));
    if p.type == :Individual
        @constraint(model, dot(p.m, s) <= p.S * p.K * (1-scale_tol));
    elseif p.type == :Total
        @constraint(model, transpose(s) * p.Q * s + dot(p.c, s) <= p.S * (1-scale_tol))
    end

    # Optimize along each sᵢ to draw a bounding box
    for sense in (MOI.MIN_SENSE, MOI.MAX_SENSE) # min and max
        for i in 1:p.K
            @objective(model, sense, s[i]);
            optimize!(model);
            if termination_status(model) == MOI.INFEASIBLE
                samp.feasible = false
                return :Infeasible
            end
            if termination_status(model) == MOI.OPTIMAL
                push!(samp.warmup, value.(s))
            end
        end
    end
    
    # Remove redundant search directions
    unique!(samp.warmup);
    samp.n_warmup = length(samp.warmup)
    if samp.n_warmup == 0
        error("Insufficient Direction Vectors Found")
    end
    samp.start = mean(samp.warmup) #? if not convex, then this start may be infeasible

    # start at the current estimated start
    samp.prev = copy(samp.start)

    samp.feasible = true
    return :Feasible
end

"""
Sample a new feasible point from the point `sampler.prev` in direction `Δ`.
"""
function step(samp::HRSampler, Δ::Vector{Float64})
    x = samp.prev
    p = samp.problem


    # No need to discuss sign of Δᵢ, all pos (neg.) λ goes to upper (lower) bounds
    # due to formal replacement of f(x) ← f(x+λΔ) and therefore given signs of coeff.

    # Apply the sᵢ>=0 constrs.
    λ_pos= -x./Δ

    # Apply the feasibility constrs.
    λ_feas = (p.N⁰ - p.Λ*x - p.c) ./ (p.Λ*Δ)

    # Apply the individual/total supply constrs.
    if p.type == :Individual
        λ_supp = (p.K*p.S - dot(p.m, x))/dot(p.m, Δ)
    elseif p.type == :Total
        A = Δ' * p.Q * Δ
        B = 2x' * p.Q * Δ + p.c' * Δ
        C = x' * p.Q * x + p.c' * x - p.S
        # @info "quad constrain with A=$(round(A, digits=3)), C=$(round(C, digits=3))"
        r0 = -B/(2A); r1 = sqrt(B^2-4A*C)/(2A)
        λ_supp = [r0-r1, r0+r1]
    end

    # ... Further extension of this work ...

    # Combine all constrs. into (λ_min to λ_max) and sample λ
    λ_all = [λ_pos; λ_supp; λ_feas]
    λ_min = isempty(λ_all[λ_all .< 0]) ? 0 : maximum(λ_all[λ_all .< 0])
    λ_max = isempty(λ_all[λ_all .> 0]) ? 0 : minimum(λ_all[λ_all .> 0])
    λ = rand(Uniform(λ_min, λ_max))

    return x + λ*Δ, λ_min, λ_max
end


function hr_sample!(samp::HRSampler, n::Int; thinning::Int = 2, burn_in::Int = -1)
    if burn_in == -1
        burn_in = round(Int, n / 2)
    end
    chain = Vector{eltype(samp.prev)}[]

    for i in 1:n
        Δ = rand(Uniform(-1, 1), samp.problem.K)
        samp.prev, λ_min, λ_max= step(samp, Δ)
        samp.n_samples += 1
        @info "Iter $i: λ sampled in range [$(round(λ_min, digits=3)), $(round(λ_max, digits=3))]"
        if i % thinning == 0 && i > burn_in
            push!(chain, copy(samp.prev))
        end
    end
    return VectorOfArray(chain)
end