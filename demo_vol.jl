
""" testing! """
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
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:Total)
end

function test_linear(K::Int=3, S::Float64=4.0, ϵ::Float64=0.75, seed::Int=42)
    Random.seed!(seed)
    σ = randn(K,K)
    c = minimum(eigvals((σ + σ')/2))
    if c < 0
        σ += (-c + ϵ) * Diagonal(fill(1,K))
    end
    Λ = inv(σ)
    Q = (Λ + Λ')/2
    m = fill(1.0,K); d = fill(1.0, K); N⁰ = fill(1e-3, K); c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:Individual)
end

function plot_sphere(sp::Sphere, itp::Union{IsoTrans, Nothing}=nothing, go_back::Bool=false)
    θ = 0:0.001:2π
    x = sp.c[1] .+ sp.r * cos.(θ)
    y = sp.c[2] .+ sp.r * sin.(θ)
    color_sp = "cyan"
    if itp !== nothing && go_back
        x = itp.invL[1,1] * x .+ itp.invL[1,2] * y
        y = itp.invL[2,1] * x .+ itp.invL[2,2] * y
        color_sp = "orange"
    end
    plot!(x, y, color=color_sp, linewidth=3)
end

function plot_samples(samples, itp::Union{IsoTrans, Nothing}=nothing, go_back::Bool=false)
    color_pts = "blue"
    if itp !== nothing && go_back
        samples = [itp.invL * y for y in samples]
        color_pts = "red"
    end
    scatter!([p[1] for p in samples], [p[2] for p in samples], color=color_pts, markersize=2, markerstrokewidth=0, ratio=1)
end


# (1) sampling linear constraint domain
p = test_linear(2, 2.5, 0.5, 32);
itp = make_isotropic(p);
region = InterPolySpheres(itp.A, itp.b, [], chevball(p, itp), :Individual);

samples = hr_sample(region, 5, 1000);
plot(ratio=1)
plot(create_poly(p), alpha = 0.25, color="purple")
plot_samples(samples, itp, true)
plot_sphere(region.chev, itp, true)


# (2) sampling quadratic constraint domain # we can lateron build wrapper function to write this region
p = test_quad(2, 20.0, 0.75);
itp = make_isotropic(p);
region = InterPolySpheres(itp.A, itp.b, [Sphere(itp.yc, sqrt(itp.t))], chevball(p, itp), :Total);
samples = hr_sample(region, 10, 2000);
plot(ratio=1)
plot_samples(samples, itp, false)
plot_sphere(region.chev, itp, false)


# (3) linear volume estimation
p = test_linear(5, 7.5, 0.5, 4);
volume_EFD(p, 10, 1, true)
volume_EFD(p, 10, 1, false)


# (4) quadratic volume estimation
p1 = test_quad(4, 8.5, 0.75, 42);
volume_EFD(p1, 10, 1, true)
volume_EFD(p1, 10, 1, false)