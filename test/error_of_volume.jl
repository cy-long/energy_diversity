"""Test for the error rate of the volume calculation, using standard examples with analytical solutions"""

include("../src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra
using ProgressMeter, Plots, SpecialFunctions

# volume of the intersection of ∑ᵢxᵢ²≤P and the positive octant of xᵢ≥0
# mimicking the omega calculation in previous results
function vol_spherical_simplex(S::Float64, K::Int)
    return (pi^(K/2) * S^(K/2)) / ((gamma(K/2 + 1) * 2^K))
end

K = 8; P_range = Vector(1.0:1.0:100.0);
vol_stan = [vol_spherical_simplex(S, K) for S in P_range]; # standard volumes
ra_stan = vol_stan[1:end-1]./vol_stan[2:end]; # standard ratios

function compute_volume(P_range, n_thread::Int, n_samples::Int, n_layers::Int)
    P_range[1] > P_range[end] || reverse!(P_range);
    
    A = -Matrix{Float64}(vcat(I(K),I(K),zeros(K)')); 
    quadbounds = [Sphere(zeros(K), sqrt(S)) for S in P_range];
    regions = [InterPolySpheres(A, vcat(zeros(2K),S), [q], :total) for (S, q) in zip(P_range, quadbounds)];
    balls = [chevball(r) for r in regions];

    v1 = volume_domain(regions[1], n_layers, 1, true);

    volumes = [0.0 for _ in P_range];
    volumes[1] = v1;
    P_range[1] > P_range[end] || reverse!(P_range)
    r = 1.0
    ratios = [];
    @showprogress desc="Volume: $(P_range[end]) to $(P_range[1])" for i in eachindex(regions)
        if r < 1e-6 || i == length(regions) || balls[i].r < 1e-9
            break
        end
        samples_i = hr_sample(regions[i], n_thread, n_samples, balls[i])
        inside_R_ii = [is_inside(s, regions[i+1]) for s in samples_i]
        r = mean(inside_R_ii)
        volumes[i+1] = volumes[i] * r
        push!(ratios, r)
    end
    reverse!(volumes); reverse!(ratios)
    P_range[end] > P_range[1] || reverse!(P_range)
    return volumes, ratios
end

# test the error in different # of samples
aa = [10, 10, 10, 10];
bb = [10^3, 10^4, 10^5, 10^6];

err_ra = Float64[]; err_va = Float64[];

for (a,b) in zip(aa,bb)
    @info "a = $a, b = $b"
    err_r = Float64[]; err_v = Float64[];
    for _ in 1:3
        volumes, ratios = compute_volume(a, b, v1, regions, balls, P_range);
        push!(err_r, norm(ra_stan - ratios)/norm(ra_stan));
        push!(err_v, norm(vol_stan - volumes)/norm(vol_stan));
    end
    push!(err_ra, mean(err_r));
    push!(err_va, mean(err_v));
end

err_ra
err_va