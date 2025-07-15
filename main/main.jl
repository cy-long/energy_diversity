""" A comprehensive survey on the feasibility domains (main results)"""
# Overview: 1. Generate a model ecosystem with S₀=8 populations with random parameters σ, d, N⁰. 2. Scale these three parameters in 3³=27 combinations; 3. downsize them by selecting [1:S] submatrix (subvector) where S∈{2,4,6,8}. 4. Compute the feasibility domains in each case. Repeat this procedure for 96 replications (in total: 96*27*4=10368 computations).

# 1. generating the model parameters: all iid. across species
# σ ∼ N(0,1²) → energetive → dissipative w/ σ₀=0.5; d ∼ N(1,0.1²); N⁰∼ N(1,0.1²)
# 2. scaling the parameters:
# σ → s_sigma * σ; s_sigma ∈{0.1, 1.0, 10.0}; d ∼ d₀ * d; d₀∈{0.1, 1.0, 10.0}; N⁰ ∼ N₀ * N⁰; N₀∈{0.1, 1.0, 10.0}
# 3. downsize the model system with S ∈ {2, 4, 6, 8} populations

include("../src/EnerFeas.jl");
using .EnerFeas
using Random, Plots, JLD2
using ProgressMeter

seed = length(ARGS) ≥ 1 ? parse(Int, ARGS[1]) : 123;
TEST_MODE = false;

# single model system test
if TEST_MODE
    p = generate_model_system(3, :total, seed, 1.0, 1.0, 1.0);
    Q_range = select_range(p);
    devols = volume_range_flux(p, Q_range);
    vols = volume_range_EFD(p, Q_range, n_sample=25000, show_p=false);
    plot(Q_range, vols./devols, xaxis=:log10,lw=1.5,xlim=(1,Q_range[end]),title="Single case")
end

# loop over 27*4 cases with 2 types
σ_scale_range = [0.1, 1.0, 10.0];
d0_range = [0.1, 1.0, 10.0];
N0_range = [0.1, 1.0, 10.0];
S_range = [2,4,6,8];
types = [:total, :indiv];
total_cases = length(σ_scale_range) * length(d0_range) * length(N0_range) * length(S_range);

if !TEST_MODE
    results = Vector{Dict}(); counter = 1;
    for (type, σsc, d0, N0) in Iterators.product(types, σ_scale_range, d0_range, N0_range)
        p0 = generate_model_system(S_range[end], :total, seed, σsc, d0, N0) # dummy, grand system
        for S in S_range
            p = sub_model_system(S, p0); Q_range = select_range(p); # inherit :total to select Q_range
            @info "Computing ($counter / $total_cases): type: $(type), S=$(S), σsc=$(σsc), d0=$(d0), N0=$(N0)\n"
            p.type = type;
            vols = volume_range_EFD(p, Q_range, n_sample=2*10^4);
            devols = volume_range_flux(p, Q_range)
            push!(results, Dict(
                :S => S, :σsc => σsc, :d0 => d0, :N0 => N0, :type => type,
                :Qs => Q_range, :vols => vols, :devols => devols,
            ))
            counter += 1;
        end
    end
end

mkpath("data/main")
@save "data/main/seed$(seed).jld2" results