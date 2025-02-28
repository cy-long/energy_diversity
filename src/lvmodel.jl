using LinearAlgebra
using Random
using Distributions

function generate_inte_matrix(Ks::Int)
    inte_matrix = -abs.(randn(Ks, Ks))
    inte_matrix[diagind(inte_matrix)] .= -1
    return inte_matrix
end


function create_sampler()
    pblm = create_energy_problem()
    Z = Matrix(Diagonal(ones(pblm.K))) # no equality const., null space=the entire space
    warmup = Vector{Vector{Float64}}()
    center = zeros(pblm.K)
    n_warmup = 0
    n_samples = 0
    prev = zeros(pblm.K)
    feas_tol = 1e-6
    bounds_tol = 1e-6
    feasible = false
    RNG = MersenneTwister(42)

    return HRSampler(pblm, Z, warmup, center, n_warmup, n_samples, prev, feas_tol, bounds_tol, feasible, RNG)
end