""" Analyze the results in theory """

using Statistics, StatsBase, Combinatorics, CategoricalArrays, Distributions
using DataFrames, Glob, JLD2
using Plots, Interpolations
using JSON

function collect_all_results(dir::String)
    files = glob("results_seed*.jld2", dir)
    all_rows = DataFrame[]
    
    for file in files
        @info "Loading $file"
        data = load(file, "results")  # Vector{Dict}

        # Extract seed from filename
        seed_match = match(r"seed(\d+)", file)
        seed = isnothing(seed_match) ? missing : parse(Int, seed_match.captures[1])

        # Add :seed key to each dictionary
        rows = [merge(d, Dict(:seed => seed)) for d in data]
        push!(all_rows, DataFrame(rows))
    end
    
    col_orders = [:seed, :S, :σsc, :d0, :N0, :type, :Qs, :vols, :devols]
    combined = vcat(all_rows...)
    return combined[!, col_orders]
end

results_df = collect_all_results("data/main")

function clean_row(row)
    return Dict(
        "S" => Int(row.S),
        "type" => String(row.type), # "total" (r) or "indiv" (i)
        "Qs" => row.Qs,
        "y"  => row.vols ./ row.devols 
    )
end

c_params = (1.0, 1.0, 1.0);
filter_cond(row) = row.σsc == c_params[1] && row.d0 == c_params[2] && row.N0 == c_params[3]
df_matr = filter(row -> filter_cond(row) && row.type == :total && row.S in [2,4,6,8], results_df)
df_init = filter(row -> filter_cond(row) && row.type == :indiv && row.S in [2,4,6,8], results_df)
data_matr = [clean_row(r) for r in eachrow(df_matr)];
data_init = [clean_row(r) for r in eachrow(df_init)];
df_single_pool = filter(row -> filter_cond(row) && row.type == :total && row.S == 2, results_df);
single_traj = isempty(df_single_pool) ? nothing : clean_row(eachrow(df_single_pool)[min(25, nrow(df_single_pool))]);

final_export = Dict(
    "data_matr" => data_matr,
    "data_init" => data_init,
    "single" => single_traj
);

open("data/output/fig3_data.json", "w") do f
    JSON.print(f, final_export)
end


df = results_r; c_params = (1.0, 1.0, 1.0); 
Gss = Vector{Vector{Float64}}();
for S in [2,4,6,8]
    i = Int(S/2)
    df_S = filter(row -> row.S==S && row.σsc==c_params[1] && row.d0==c_params[2] && row.N0==c_params[3], df);
    Gs = Vector{Float64}();
    for row in eachrow(df_S)
        push!(Gs, score_unimodal(row.Qs, row.vols ./ row.devols))
    end
    push!(Gss, Gs)
end

mean(Gss[1])
mean(Gss[2])
mean(Gss[3])
mean(Gss[4])