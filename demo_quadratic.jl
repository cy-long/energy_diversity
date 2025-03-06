using Revise

using Random
using LinearAlgebra
using Distributions
using Plots

includet("src/lvmodel.jl")
includet("src/sampler.jl")

const GRB_ENV = Gurobi.Env(output_flag=0)

p = create_energy_problem(4, 3.5, :Total, seed=467);
make_problem_dissipative!(p)
samp = create_sampler(p);

# test: warmup works for quadratic cases
warmup!(samp)


# test for sampling
samples = Vector{Vector{Float64}}()

for i in 1:5000
    Δ = rand(Uniform(-1, 1), p.K)
    samp.prev, λmin, λmax= step(samp, Δ)
    samp.n_samples += 1
    push!(samples, copy(samp.prev))
end

supply = [total_supply(sample, p) for sample in samples]

supply = [individual_supply(sample, p) for sample in samples]


plot(1:5000, supply, xlabel="sampling run", ylabel="sampling supply", label="S")
