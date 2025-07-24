""" A comprehensive survey on the feasibility domains (main results)"""
# Overview: 1. Generate a model ecosystem with S₀=8 populations with random parameters σ, d, N⁰. 2. Scale these three parameters in 9 combinations; 3. downsize them by selecting [1:S] submatrix (subvector) where S ∈ {2,4,6,8}. 4. Compute the feasibility domains in each case. Repeat this procedure for 100 replications.

# 1. generating the model parameters: all iid. across species
# σ ∼ N(0,1²) → energetive → dissipative w/ σ₀=0.5; d ∼ N(1,0.1²); N⁰∼ N(1,0.1²)
# 2. scaling the parameters:
# σ → s_sigma * σ; s_sigma ∈{0.5, 1.0, 2.0}; d ∼ d₀ * d; d₀∈{0.5, 1.0, 2.0}; N⁰ ∼ N₀ * N⁰; N₀∈{0.5, 1.0, 2.0}
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
    vols = volume_range_EFD(p, Q_range, n_sample=2500, show_p=true);
    plot(Q_range, vols./devols, xaxis=:log10,lw=1.5,xlim=(1,Q_range[end]),title="Single case")
    savefig("figure/test_single_case.pdf")
end

@info "Running main computation with seed=$(seed)"
params_range = [
    (1.0, 1.0, 1.0), 
    (1.0, 2.0, 1.0), (1.0, 0.5, 1.0),
    (0.5, 1.0, 1.0), (2.0, 1.0, 1.0),
    (1.0, 1.0, 0.5), (1.0, 1.0, 2.0),
    (2.0, 2.0, 1.0), (0.5, 0.5, 1.0)
];
S_range = [2, 4, 6, 8];
types = [:total, :indiv];
total_cases = length(params_range) * length(S_range) * length(types);

function run_main(seed)
    results = Vector{Dict}(); counter = 1;
    for (type, (σsc, d0, N0)) in Iterators.product(types, params_range)
        p0 = generate_model_system(S_range[end], :total, seed, σsc, d0, N0) # dummy, grand system
        for S in S_range
            p = sub_model_system(S, p0); Q_range = select_range(p);
            @info "Computing ($counter / $total_cases): type: $(type), S=$(S), σsc=$(σsc), d0=$(d0), N0=$(N0)\n"
            p.type = type;
            vols = volume_range_EFD(p, Q_range, n_sample=2*10^4, show_p=true, show_dt = 12.0);
            devols = volume_range_flux(p, Q_range)
            push!(results, Dict(
                :S => S, :σsc => σsc, :d0 => d0, :N0 => N0, :type => type,
                :Qs => Q_range, :vols => vols, :devols => devols,
            ))
            counter += 1;
            if TEST_MODE && counter > 2
                @info "Test mode: stopping after 10 cases"
                return results;
            end
        end
    end
    return results
end

results = run_main(seed);

mkpath("data/main")
@save "data/main/results_seed$(seed).jld2" results