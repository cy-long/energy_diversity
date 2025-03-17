"""
This script develops tools to sample on energetic feasibility domains. It is mainly based
on hit-and-run sampling on convex set. Which is inspired by COBRA¹ and Biologically Constrained Foodwebs²

¹ https://opencobra.github.io/cobrapy
² https://doi.org/10.1073/pnas.2212061120
"""

using LinearAlgebra
using Random
using Distributions
using RecursiveArrayTools
using JuMP
using Gurobi

include("_problem.jl")


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


function create_model(p::EnergyConstraintProblem, abs_tol::Float64=0.0, scale_tol::Float64=0.0)
    model = Model()
    @variable(model, s[i = 1:p.K] >= abs_tol);
    @constraint(model, p.Λ * s .>= (p.N⁰ - p.c) * (1 + scale_tol));
    flag_lin = false
    if p.type == :Individual
        @constraint(model, dot(p.m, s) <= p.S * p.K * (1-scale_tol));
        flag_lin = true
    elseif p.type == :Total
        @constraint(model, transpose(s) * p.Q * s + dot(p.c, s) <= p.S * (1-scale_tol))
    end
    return model, s, flag_lin
end


function warmup!(samp::HRSampler)
    p = samp.problem
    abs_tol= samp.abs_tol
    scale_tol = samp.scale_tol
    
    # Initialize a JuMP Model with Gurobi solver for the warmup
    model, s, _ = create_model(p, abs_tol, scale_tol)
    set_optimizer(model, () -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(model, "OutputFlag", 0)   # Mute all Gurobi output
    set_optimizer_attribute(model, "LogToConsole", 0) # Prevent logging to console

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
function hr_step(samp::HRSampler, Δ::Vector{Float64})
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


function hr_sample!(samp::HRSampler, n::Int; thinning::Int = 1, burn_in::Int = -1)
    if burn_in == -1
        burn_in = round(Int, n / 2)
    end

    chain = []
    for i in 1:n
        Δ = rand(Uniform(-1, 1), samp.problem.K)
        samp.prev, λ_min, λ_max= hr_step(samp, Δ)
        samp.n_samples += 1
        if i % 500 == 0
            @info "Iter $i: λ sampled in range [$(round(λ_min, digits=3)), $(round(λ_max, digits=3))]"
        end
        if i % thinning == 0 && i > burn_in
            push!(chain, copy(samp.prev))
        end
    end
    return VectorOfArray(chain)
end