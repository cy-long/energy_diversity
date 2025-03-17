"""
This script provides tools to sample and measure quadratic EFD. Methods include Cholesky decomposition of the quadratic boudnary (eclipsoid) and finding Chebyshev ball to warm up sampling and set reference to the volume estimation.

Ref:
¹https://doi.org/10.1016/j.comgeo.2022.101916
²https://doi.org/10.1145/3194656
"""
# some parts should be extracted into sampler.jl. e.g. chevball warmup

using Random, Distributions
using JuMP, Gurobi
using Plots
using LinearAlgebra
using SpecialFunctions
using .EnergyProblem

struct Sphere
    c::Vector{Float64} # center of the sphere
    r::Float64 # radius of the sphere
end

function vol_sphere(sp::Sphere)
    K = length(sp.c)
    return (pi^(K/2) / gamma(K/2 + 1)) * sp.r^K
end

function isotropic_transform(p::EnergyConstraintProblem)
    p.type == :Total || error("Only Total problems are supported")
    L = cholesky(p.Q).U
    A = -vcat(I, p.Λ) * inv(L)
    b = -vcat(zeros(p.K), p.N⁰ - p.c)
    yc = -0.5 * inv(L)' * p.c
    t = p.S + 0.25 * p.c' * inv(p.Q) * p.c
    invL = inv(L)
    return IsotropicTransformParams(L, invL, yc, t, A, b)
end

# which will replace the warm-up stage
#? add feasible indicator
function chevball(p::EnergyConstraintProblem, itp::IsotropicTransformParams)
    yc = itp.yc; A = itp.A; b = itp.b; t = itp.t
    A_norms = [norm(row) for row in eachrow(A)]
    model = Model()
    @variable(model, x[i = 1:p.K])
    @variable(model, r >= 0)
    @constraint(model, A * x + A_norms * r .<= b)
    @constraint(model, [sqrt(t)-r; x - yc] in SecondOrderCone())
    @objective(model, Max, r)
    set_optimizer(model, () -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "LogToConsole", 0)
    optimize!(model)
    return Sphere(value.(x), value.(r))
end

function hr_step_quad(x::Vector{Float64}, itp::IsotropicTransformParams, sp::Union{Sphere, Nothing}=nothing)
    yc = itp.yc; A = itp.A; b = itp.b; t = itp.t
    Δ = randn(size(yc,1))
    λ_lin = (b - A*x) ./ (A * Δ)

    a = dot(Δ, Δ); b = 2*dot(Δ, x-yc); c = dot(x-yc,x-yc) - t; d = sqrt(b^2 - 4*a*c)
    λ_quad = [(-b-d)/(2*a), (-b+d)/(2*a)]

    λ_all = [λ_lin; λ_quad]

    if !isnothing(sp)
        a1 = dot(Δ, Δ); b1 = 2*dot(Δ, x-sp.c); c1 = dot(x-sp.c,x-sp.c)-sp.r^2; d1 = sqrt(b1^2 - 4*a1*c1)
        λ_quad1 = [(-b1-d1)/(2*a1), (-b1+d1)/(2*a1)]
        λ_all = vcat(λ_all, λ_quad1)
    end

    λ_min = isempty(λ_all[λ_all .< 0]) ? 0 : maximum(λ_all[λ_all .< 0])
    λ_max = isempty(λ_all[λ_all .> 0]) ? 0 : minimum(λ_all[λ_all .> 0])
    λ = rand(Uniform(λ_min, λ_max))

    return x + λ * Δ, λ_min, λ_max
end


# reconstuct the hr sampler function, since we can use the chevball as warmup
function hr_sample_quad(
    p::EnergyConstraintProblem,
    n_threads::Int=4,
    n_samples::Int=2500,
    go_back::Bool=false,
    itp::IsotropicTransformParams=isotropic_transform(p);
    chev::Sphere=chevball(p, itp),
    sp::Union{Sphere, Nothing} = nothing
    )

    # itp = isotropic_transform(p) # could be called as a parameter
    samples = []
    
    for _ in 1:n_threads
        e = randn(p.K)
        x = chev.c + chev.r * rand(Uniform(0,1)) * e / norm(e)
        thread = []
        for _ in 1:n_samples
            x, _, _ = hr_step_quad(x, itp, sp)
            push!(thread, x)
        end
    samples = vcat(samples, thread)
    end
    
    if go_back
        samples = [itp.invL * y for y in samples]
    end
    return(samples)
end

function isin_domain_and_sphere(
    x::Vector{Float64},
    itp::IsotropicTransformParams,
    sp::Sphere)
    return all(itp.A * x .<= itp.b) && norm(x - itp.yc) <= sqrt(itp.t) && norm(x - sp.c) <= sp.r
end


#! seems problematic now. threads implementation, N/eN hyper parameters, skipping W steps
function volume_quad(p::EnergyConstraintProblem, itp::IsotropicTransformParams, go_back::Bool=false, N::Int=10, eN::Int=1)
    chev_p = chevball(p, itp)
    
    # briefly explore the domain to estimate ρ
    samples_0 = hr_sample_quad(p, 1, 5000, false, itp, chev=chev_p);
    r = chev_p.r
    ρ = maximum([norm(x-chev_p.c) for x in samples_0])

    # Multiphase MC
    rads = [r*(ρ/r)^(k/N) for k in 0:(N+eN)]
    vol_ratio = Float64[]
    for i in eachindex(rads)
        if i == firstindex(rads)
            push!(vol_ratio, vol_sphere(Sphere(chev_p.c, rads[i])))
            continue
        end
        sphere_focu = Sphere(chev_p.c, rads[i])
        sphere_prev = Sphere(chev_p.c, rads[i-1])
        # plot()
        # plot_sphere(sphere_focu)
        # plot_sphere(sphere_prev)
        samples_i = hr_sample_quad(p, 10, 1000, false, itp, chev=chev_p, sp=sphere_focu)
        # scatter!([s[1] for s in samples_i], [s[2] for s in samples_i], color="cyan", markersize=1, markerstrokewidth=0, ratio=1)
        push!(vol_ratio, mean([isin_domain_and_sphere(x, itp, sphere_prev) for x in samples_i]))
    end
    if go_back
        push!(vol_ratio, det(itp.invL))
    end
    return(prod(vol_ratio))
end


# ---- only for temporarily testing ----
function test_quad(K::Int=3, S::Float64=4.0, ϵ::Float64=0.75, seed::Int=42)
    Random.seed!(seed)
    σ = randn(K,K)
    # σ = Matrix(1.0I, K, K)
    c = minimum(eigvals((σ + σ')/2))
    if c < 0
        σ += (-c + ϵ) * Diagonal(fill(1,K))
    end
    Λ = inv(σ)
    Q = (Λ + Λ')/2
    m = fill(1.0,K); d = fill(1.0, K); N⁰ = fill(1e-3, K); c = -Λ * d
    return EnergyConstraintProblem(σ,Λ,Q,c,m,d,N⁰,K,S,:Total)
end

include("EnergyProblem.jl");
const GRB_ENV = Gurobi.Env(output_flag=0);

p = test_quad(6, 16.0, 0.5, 42);
itp = isotropic_transform(p);
chev_p = chevball(p, itp);

# Chebyshev warmup + sampling in y-space
# samples = hr_sample_quad(p, 4, 2500)

# pt_test = rand(samples)
# isin_domain_and_sphere(pt_test, itp, xc_chev, 2*r_chev)

function plot_sphere(sp::Sphere)
    θ = 0:0.001:2π
    x = sp.c[1] .+ sp.r * cos.(θ)
    y = sp.c[2] .+ sp.r * sin.(θ)
    plot!(x, y, color="cyan", linewidth=1.5)
end


volume_quad(p, itp, 10, 1)

volume_quad(p, itp, 20, 1)

volume_quad(p, itp, 20, 2)

volume_quad(p, itp, 5, 2)


boundary_s = sqrt(itp.t)* [[cos(θ), sin(θ)] for θ in 0:0.001:2π] .|> 
    (x -> itp.invL * x + itp.invL * itp.yc) |>
    x -> sort(x, by = s -> s[1], rev = true) |>
    x -> filter(s -> s[1] >= 0 && s[2] >= 0, x);
chev_x = [xc_chev .+ r_chev .* [cos(θ), sin(θ)] for θ in 0:0.001:2π];
chev_s = [itp.invL * x for x in chev_x];


# Plotting in s-space
scatter([s[1] for s in samples], [s[2] for s in samples], color="blue", markersize=1, markerstrokewidth=0, label="Samples")
plot!([s[1] for s in boundary_s], [s[2] for s in boundary_s], color="blue", linewidth=2, label="Tot. Supp. Bound")
plot!([s[1] for s in chev_s], [s[2] for s in chev_s], color="red", linewidth=3, label="Chebyshev Ball")

# Plotting in y-space
scatter([((itp.L)*s)[1] for s in samples], [((itp.L)*s)[2] for s in samples], color="darkorange", markersize=1, markerstrokewidth=0, label="Samples in y", ratio=1)
plot!([y[1] for y in chev_x], [y[2] for y in chev_x], color="blue", linewidth=2, label="Chebyshev Ball")