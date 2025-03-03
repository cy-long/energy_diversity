using LinearAlgebra
using Random
using Distributions
using RecursiveArrayTools

using JuMP
using Gurobi
ENV["GRB_LICENSE_FILE"] = "/Users/longcy/Documents/Gurobi/gurobi.lic" 
# replace with your Gurobi license path
const GRB_ENV = Gurobi.Env(output_flag=0)

include("lvmodel.jl")


function warmup_sampler!(samp::HRSampler)
    b_tol= samp.bounds_tol
    f_tol = samp.feas_tol
    N_eff = samp.pblm.Λ * samp.pblm.d + samp.pblm.N⁰
    
    # Initialize a JuMP Model with Gurobi as the solver
    model = Model(() -> Gurobi.Optimizer(GRB_ENV));
    set_optimizer_attribute(model, "OutputFlag", 0)   # Mutes all Gurobi output
    set_optimizer_attribute(model, "LogToConsole", 0) # Prevents logging to console
    # set_optimizer_attribute(model, "Threads", 1)

    # setup the warmup problem
    @variable(model, s[i = 1:samp.pblm.K] >= b_tol);
    @constraint(model, samp.pblm.Λ * s .>= N_eff * (1+f_tol));
    @constraint(model, dot(samp.pblm.m, s) <= samp.pblm.Sᵢ * samp.pblm.K* (1-b_tol)) ;

    # Warm up along each individual dimension to get the Min and Max to form a spanning set
    for sense in (MOI.MIN_SENSE, MOI.MAX_SENSE)
        for i in 1:samp.pblm.K
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
    samp.n_samples = samp.n_warmup
    samp.center = mean(samp.warmup) # if not convex, then this center may be infeasible

    # start at the current estimated center
    samp.prev = copy(samp.center)

    samp.feasible = true
    return :Feasible
end

"""
Sample a new feasible point from the point `sampler.prev` in direction `Δ`.
"""
function hr_step(samp::HRSampler, Δ::Vector{Float64})
    x = samp.prev
    b_tol = samp.bounds_tol
    f_tol = samp.feas_tol
    Λ = samp.pblm.Λ
    N⁰ = samp.pblm.N⁰
    m = samp.pblm.m
    K = samp.pblm.K
    d = samp.pblm.d
    Sᵢ = samp.pblm.Sᵢ
    # dir=1 if Δ>feas_tol; dir=-1 if Δ<-feas_tol; dir=0 elsewise
    direction = sign.(Δ) .* (abs.(Δ) .>= f_tol)
    valid = findall(direction .!= 0)
    
    if valid == []
        pass # I don't know what to do here
    end

    # sᵢ>=0
    λ1= -(1-b_tol) * x[valid]./Δ[valid] # positive goes to upper, negative goes to lower

    #∑ᵢ sᵢmᵢ <= K*Sᵢ
    λ2 = (K*Sᵢ - dot(x, m))/dot(Δ, m) #positive goes to upper, negative goes to lower

    # Λ(s-d) >= N⁰
    λ3 = (N⁰ - Λ*x + Λ*d) ./ (Λ*Δ) # positive goes to upper, negative goes to lower

    # combine 1,2,3, find the miminal of the positive and the maximum of the negative
    λall = [λ1; λ2; λ3]
    λmin = maximum(λall[λall .< 0])
    λmax = minimum(λall[λall .> 0])
    # sample λ from (λmin to λmax)
    λ = rand(Uniform(λmin, λmax))
    return x + λ * Δ
end


function hr_sample!(samp::HRSampler, n::Int; thinning::Int = 1, burn_in::Int = -1)
    if burn_in == -1 # what is the correct way to do this in a typed manner?
        burn_in = round(Int, n / 2)
    end
    chain = Vector{eltype(samp.prev)}[]

    for i in 1:n
        Δ = samp.Z * rand(Uniform(-1, 1), samp.pblm.K)
        samp.prev = hr_step(samp, Δ)
        samp.n_samples += 1
        if i % thinning == 0 && i > burn_in
            push!(chain, copy(samp.prev))
        end
    end
    return VectorOfArray(chain);

end