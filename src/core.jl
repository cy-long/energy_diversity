"""
This scripts implements the backend computational geometry. Constrained energetic problem can be  generalized into the InterPolySpheres type, where
- The convex polyhedron is from feasibility conditions;
- The spheres are total energetic boundary (if applicable) and auxillary spheres.

Around InterPolySpheres we build sampling and volume estimation methods. Sampling is perfomed by by RDHR (random direction hit-and-run¹²), Volume estimation is performed by MMC (multiphase Monte Carlo³⁴).

For linear constraints, trianglation method is included as a comparison.

References: ¹https://opencobra.github.io/cobrapy, ²https://doi.org/10.1073/pnas.2212061120
³https://doi.org/10.1016/j.comgeo.2022.101916, ⁴https://doi.org/10.1145/3194656
"""
function min_vol_ellipsoid(po::Polyhedron)
    pts = hcat(points(po)...)
    ϵ = minimum_volume_ellipsoid(pts)
    return Matrix(ϵ.H)
end

struct Sphere
    c::Vector{Float64} # center
    r::Float64 # radius
end
function vol_sphere(sp::Sphere)
    K = length(sp.c)
    return (pi^(K/2) / gamma(K/2 + 1)) * sp.r^K
end

struct InterPolySpheres
    A::Matrix{Float64}
    b::Vector{Float64}
    sps::Vector{Sphere}
    type::Symbol # :Individual or :Total
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
        error("Chevball optimization failed: $(termination_status(model))")
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
    λ = rand(Uniform(λ_min, λ_max))

    return x + λ * Δ, λ_min, λ_max
end

# regions will be intersecting domain with k*chevball(domain), whose chevball is the same as for the domain
# therefore, no need to recompute chevball for each region, as it could be costy in time
function hr_sample(region::InterPolySpheres, n_threads::Int=10, n_samples::Int=2000, chev_start::Sphere=chevball(region))
    K = length(chev_start.c)
    samples = []
    burn_in = floor(Int, n_samples/2) # burn-in period

    for _ in 1:n_threads
        e = randn(K)
        # starting from random point inside the chevball
        x = chev_start.c + rand(Uniform(0, chev_start.r/norm(e))) * e
        thread = []
        for i in 1:n_samples
            x, _, _ = hr_step(x, region)
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

    # create and sandwhich the domain
    if domain.type == :Individual
        po = polyhedron(hrep(domain.A, domain.b))
        if exact
            return Polyhedra.volume(po)
        end
        vertices = points(po)
        ρ = maximum([norm(v-chev.c) for v in vertices])
    elseif domain.type == :Total
        sp0 = domain.sps[1]
        samples_0 = hr_sample(domain, 1, 5000, chev)
        ρ = maximum([norm(x-chev.c) for x in samples_0]) # slightly dumb way to sandwhich the domain
    end

    # slice the domain with eccentric spheres
    r_phs = [r*(ρ/r)^(k/N) for k in 0:(N+eN)]
    sp_phs = [Sphere(chev.c, r_ph) for r_ph in r_phs]

    if domain.type == :Individual 
        regions = [InterPolySpheres(domain.A, domain.b, [sp], :Individual) for sp in sp_phs]
    elseif domain.type == :Total
        regions = [InterPolySpheres(domain.A, domain.b, [sp0, sp], :Total) for sp in sp_phs]
    end

    # perform multiphase Monte Carlo sampling
    vol_ratio = Float64[]
    for i in eachindex(r_phs)
        if i == firstindex(r_phs)
            push!(vol_ratio, vol_sphere(sp_phs[i]))
            continue
        end
        samples_i = hr_sample(regions[i], 10, 5000, chev) #TODO: start chains from previous samples, instead of chevball
        inside_i = [is_inside(x, regions[i-1]) for x in samples_i]
        push!(vol_ratio, 1 / mean(inside_i))
    end

    return prod(vol_ratio)
end