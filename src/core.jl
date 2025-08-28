"""
This scripts implements the backend computational geometry. Constrained energetic problem can be  generalized into the InterPolyBalls type, where
- The convex polyhedron is from feasibility conditions;
- The spheres are total energetic boundary (if applicable) and auxillary spheres.

Around InterPolyBalls we build sampling and volume estimation methods. Sampling is perfomed by by RDHR (random direction hit-and-run¹²), Volume estimation is performed by MMC (multiphase Monte Carlo³⁴).

For linear constraints, trianglation method is included as a comparison.

References: ¹https://opencobra.github.io/cobrapy, ²https://doi.org/10.1073/pnas.2212061120
³https://doi.org/10.1016/j.comgeo.2022.101916, ⁴https://doi.org/10.1145/3194656
"""

struct Sphere
    c::Vector{Float64} # center
    r::Float64 # radius
end

struct InterPolyBalls
    A::Matrix{Float64}
    b::Vector{Float64}
    sps::Vector{Sphere}
    type::Symbol # :indiv or :total
end

function rand_sphere(sp::Sphere, cut::Float64)
    0.0 < cut <= 1.0 || throw(ArgumentError("cut must be in (0, 1]"))
    K = length(sp.c)
    e = randn(K)
    x = sp.c + rand(Uniform(0, cut * sp.r/norm(e))) * e
    return x
end

function rand_warmup(domain::InterPolyBalls, ball::Sphere, n_threads::Int, max_trials::Int)
    x_warmup = Vector{Vector{Float64}}(undef, n_threads)
    k = 0
    for _ in 1:max_trials
        x = rand_sphere(ball, 0.95)
        if is_inside(x, domain)
            x_warmup[k+1] = x
            k += 1
            if k >= n_threads
                return x_warmup
            end
        end
    end
    error("rand_warmup: Failed to generate $n_threads seeds after $max_trials trials.")
end

function vol_sphere(sp::Sphere)
    K = length(sp.c)
    return (pi^(K/2) / gamma(K/2 + 1)) * sp.r^K
end

function is_inside_sphere(x::Vector{Float64}, sp::Sphere)
    return norm(x - sp.c) <= sp.r
end

function chevball(domain::InterPolyBalls)
    A = domain.A; b = domain.b; K=size(A,2)
    A_norms = [norm(row) for row in eachrow(A)]

    model = Model()
    @variable(model, x[i = 1:K])
    @variable(model, r >= 0)
    @constraint(model, A * x + A_norms * r .<= b)
    for sp in domain.sps
        @constraint(model, [sp.r - r; x - sp.c] in SecondOrderCone())
    end
    @objective(model, Max, r)
    
    set_optimizer(model, SCS.Optimizer); #> handling usage outside this package
    set_optimizer_attribute(model, "verbose", 0)
    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        return Sphere(fill(NaN, K), 0.0) # indicate infeasible but not crash
    else
        return Sphere(value.(x), value.(r))
    end
end

function hr_step(x::Vector{Float64}, region::InterPolyBalls)
    Δ = randn(size(x,1))
    A = region.A; b = region.b
    λs = (b - A*x) ./ (A * Δ)

    for sp in region.sps
        yc = sp.c; r = sp.r
        a = dot(Δ, Δ); b = 2*dot(Δ, x-yc); c = dot(x-yc,x-yc)-r^2; 
        if b^2 - 4*a*c <= 0
            @info "HR stuck with ball boundaries"
            return nothing
        end
        d = sqrt(b^2 - 4*a*c)
        λs = vcat(λs, [(-b-d)/(2*a), (-b+d)/(2*a)])
    end

    λ_min = isempty(λs[λs .< 0]) ? 0 : maximum(λs[λs .< 0])
    λ_max = isempty(λs[λs .> 0]) ? 0 : minimum(λs[λs .> 0])
    if λ_min < λ_max-1e-6
        λ = rand(Uniform(λ_min, λ_max))
        return x + λ * Δ
    else
        @info "HR stuck with polyhedron boundaries"
        return nothing
    end
end

# regions will be intersecting domain with ϵ*chevball(domain), whose chevball is the same as for the domain
# therefore, no need to recompute chevball for each region, as it could be costy in time
function hr_sample(region::InterPolyBalls, n_threads::Int=10, n_samples::Int=2000, start::Union{Sphere,Vector{Vector{Float64}}}=chevball(region))
    if isa(start, Sphere)
        x_seeds = rand_warmup(region, start, n_threads, 100*n_threads)
    else
        size(start,1) >= n_threads || throw(ErrorException("insufficient sampling starting points"))
        x_seeds = start
    end

    burn_in = floor(Int, n_samples/2)
    idx = 1; samples = Vector{Vector{Float64}}(undef, n_threads * (n_samples - burn_in))

    for i in 1:n_threads # can be parallelized
        x = x_seeds[i]
        for j in 1:n_samples
            x1 = hr_step(x, region)
            if isnothing(x1)
                if j <= burn_in
                    @warn "HR sampling gets stuck in early steps"
                end
                x1 = hr_step(x_seeds[i], region) # manually reset to start
            end
            if j > burn_in
                samples[idx] = x1; idx += 1; x = x1
            end
        end
    end
    return samples
end

function is_inside(x::Vector{Float64}, region::InterPolyBalls)
    in_sps = [norm(x - sp.c) <= sp.r for sp in region.sps]
    in_po = region.A * x .<= region.b
    return all(in_sps) && all(in_po)
end

function volume_domain(domain::InterPolyBalls, N::Int=10, eN::Int=1, exact::Bool=false)
    chev = chevball(domain)
    r = chev.r
    if r == 0.0 
        @warn "domain has zero volume: empty chevball"
        return 0.0
    end

    # create and sandwhich the domain
    if domain.type == :indiv
        po = polyhedron(hrep(domain.A, domain.b))
        if exact
            return Polyhedra.volume(po)
        end
        vertices = points(po)
        ρ = maximum([norm(v-chev.c) for v in vertices])
    elseif domain.type == :total
        sp0 = domain.sps[1]
        samples_0 = hr_sample(domain, 1, 10000, chev)
        ρ = maximum([norm(x-chev.c) for x in samples_0]) # slightly dumb way to sandwhich the domain
    end

    # slice the domain with eccentric spheres
    r_phs = [r*(ρ/r)^(k/N) for k in 0:(N+eN)]
    sp_phs = [Sphere(chev.c, r_ph) for r_ph in r_phs]

    if domain.type == :indiv 
        regions = [InterPolyBalls(domain.A, domain.b, [sp], :indiv) for sp in sp_phs]
    elseif domain.type == :total
        regions = [InterPolyBalls(domain.A, domain.b, [sp0, sp], :total) for sp in sp_phs]
    end

    # perform multiphase Monte Carlo sampling
    vol_ratio = Vector{Float64}(undef, length(r_phs))
    for i in eachindex(r_phs)
        if i == firstindex(r_phs)
            vol_ratio[i] = vol_sphere(sp_phs[i])
            continue
        end
        samples_i = hr_sample(regions[i], 10, 10000, chev)
        inside_i = [is_inside(x, regions[i-1]) for x in samples_i]
        vol_ratio[i] = 1 / mean(inside_i)
    end

    return prod(vol_ratio)
end