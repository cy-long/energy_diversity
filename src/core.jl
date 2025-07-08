"""
This scripts implements the backend computational geometry. Constrained energetic problem can be  generalized into the InterPolySpheres type, where
- The convex polyhedron is from feasibility conditions;
- The spheres are total energetic boundary (if applicable) and auxillary spheres.

Around InterPolySpheres we build sampling and volume estimation methods. Sampling is perfomed by by RDHR (random direction hit-and-run¹²), Volume estimation is performed by MMC (multiphase Monte Carlo³⁴).

For linear constraints, trianglation method is included as a comparison.

References: ¹https://opencobra.github.io/cobrapy, ²https://doi.org/10.1073/pnas.2212061120
³https://doi.org/10.1016/j.comgeo.2022.101916, ⁴https://doi.org/10.1145/3194656
"""

struct Sphere
    c::Vector{Float64} # center
    r::Float64 # radius
end

function rand_sphere(sp::Sphere)
    K = length(sp.c)
    e = randn(K)
    x = sp.c + rand(Uniform(0, sp.r/norm(e))) * e
    return x
end

function vol_sphere(sp::Sphere)
    K = length(sp.c)
    return (pi^(K/2) / gamma(K/2 + 1)) * sp.r^K
end

function is_inside_sphere(x::Vector{Float64}, sp::Sphere)
    return norm(x - sp.c) <= sp.r
end

struct InterPolySpheres
    A::Matrix{Float64}
    b::Vector{Float64}
    sps::Vector{Sphere}
    type::Symbol # :indiv or :total
end

function chevball(domain::InterPolySpheres)
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
    
    set_optimizer(model, () -> Gurobi.Optimizer(GRB_ENV)); #> handling usage outside this package
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "LogToConsole", 0)
    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        return Sphere(fill(NaN, K), 0.0) # indicate infeasible but not crash
    else
        return Sphere(value.(x), value.(r))
    end
end

function hr_step(x::Vector{Float64}, region::InterPolySpheres)
    Δ = randn(size(x,1))

    A = region.A; b = region.b
    λs = (b - A*x) ./ (A * Δ)

    for sp in region.sps
        yc = sp.c; r = sp.r
        a = dot(Δ, Δ); b = 2*dot(Δ, x-yc); c = dot(x-yc,x-yc)-r^2; d = sqrt(b^2 - 4*a*c)
        λs = vcat(λs, [(-b-d)/(2*a), (-b+d)/(2*a)])
    end

    λ_min = isempty(λs[λs .< 0]) ? 0 : maximum(λs[λs .< 0])
    λ_max = isempty(λs[λs .> 0]) ? 0 : minimum(λs[λs .> 0])
    if λ_min < λ_max-1e-7
        λ = rand(Uniform(λ_min, λ_max))
        return x + λ * Δ
    else
        return nothing
    end
end

# regions will be intersecting domain with ϵ*chevball(domain), whose chevball is the same as for the domain
# therefore, no need to recompute chevball for each region, as it could be costy in time
function hr_sample(region::InterPolySpheres, n_threads::Int=10, n_samples::Int=2000, 
    start::Union{Sphere,Vector{Vector{Float64}}}=chevball(region))
    if isa(start, Sphere)
        seeds = [rand_sphere(start) for _ in 1:n_threads]    
    else
        size(start,1) >= n_threads || throw(ErrorException("insufficient sampling seeds"))
        seeds = start
    end
    
    samples = []
    burn_in = floor(Int, n_samples/2) # burn-in period

    for i in 1:n_threads
        x = seeds[i]
        thread = []
        for i in 1:n_samples
            x = hr_step(x, region)
            if isnothing(x)
                x = hr_step(seeds[i], region)
                @warn "HR failed, restart at seeds"
            end
            if i > burn_in
                push!(thread, x) #TODO: thinning?
            end
        end
    samples = vcat(samples, thread)
    end
    return samples
end

function is_inside(x::Vector{Float64}, region::InterPolySpheres)
    in_sps = [norm(x - sp.c) <= sp.r for sp in region.sps]
    in_po = region.A * x .<= region.b
    return all(in_sps) && all(in_po)
end

function volume_domain(domain::InterPolySpheres, N::Int=10, eN::Int=1, exact::Bool=false)
    chev = chevball(domain)
    r = chev.r
    if r == 0.0 
        @warn "Chevball is infeasible in volume_domain"
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
        samples_0 = hr_sample(domain, 1, 5000, chev)
        ρ = maximum([norm(x-chev.c) for x in samples_0]) # slightly dumb way to sandwhich the domain
    end

    # slice the domain with eccentric spheres
    r_phs = [r*(ρ/r)^(k/N) for k in 0:(N+eN)]
    sp_phs = [Sphere(chev.c, r_ph) for r_ph in r_phs]

    if domain.type == :indiv 
        regions = [InterPolySpheres(domain.A, domain.b, [sp], :indiv) for sp in sp_phs]
    elseif domain.type == :total
        regions = [InterPolySpheres(domain.A, domain.b, [sp0, sp], :total) for sp in sp_phs]
    end

    # perform multiphase Monte Carlo sampling
    vol_ratio = Float64[]
    for i in eachindex(r_phs)
        if i == firstindex(r_phs)
            push!(vol_ratio, vol_sphere(sp_phs[i]))
            continue
        end
        samples_i = hr_sample(regions[i], 20, 10000, chev)
        inside_i = [is_inside(x, regions[i-1]) for x in samples_i]
        push!(vol_ratio, 1 / mean(inside_i))
    end

    return prod(vol_ratio)
end