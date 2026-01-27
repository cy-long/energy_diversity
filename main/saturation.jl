using EnerFeas
using LinearAlgebra, Random, Distributions
using DataFrames
using JLD2, JSON

S = 2;
Qs = 10.0 .^ range(0, stop=3, length=300);
ks = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5]

vols_all = [];
for seed in 1:25
    @info "Computing for seed $seed"
    p = generate_model_system(S, seed, 1.0, 1.0, 1.0);
    vols_c = volume_range_C(p, Qs);
    vols_seed = [];
    for kval in ks;
        p.k = kval * ones(S);
        vols = volume_range_EFD(p, :matr, Qs);
        push!(vols_seed, vols);
    end
    push!(vols_all, vols_seed);
end

rows_buffer = []
for (seed, vols_seed) in enumerate(vols_all)
    p = generate_model_system(S, seed, 1.0, 1.0, 1.0)
    vols_c = volume_range_C(p, Qs)
    for (i, vols) in enumerate(vols_seed)
        push!(rows_buffer, (seed = seed, k = ks[i], Qs = Qs, vols = vols, vols_c = vols_c))
    end
end

df_all = DataFrame(rows_buffer)
@save "result_saturation-S$(S).jld2" vols_all Qs ks df_all

S = 2;
@load "result_saturation-S$(S).jld2"
export_data = [
    Dict(
        "k" => row.k,
        "Qs" => row.Qs,
        "y" => row.vols ./ row.vols_c
    ) 
    for row in eachrow(df_all)
];

open("data/output/k_curves_data-S$(S).json", "w") do f
    JSON.print(f, export_data)
end