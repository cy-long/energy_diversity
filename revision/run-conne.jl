""" Sensitivity of connectance on PM """

using EnerFeas
using JLD2

# ARGS: seed [sampling_chains] [n_sample] [n_layer] [show_prog] [outdir]

function parse_arg(args::Vector{String}, idx::Int, default, ::Type{T}) where {T}
    return length(args) >= idx ? parse(T, args[idx]) : default
end

function unique_outfile(outdir::AbstractString, seed::Int)
    base = joinpath(outdir, "connectance_seed$(seed)")
    outfile = base * ".jld2"
    k = 1
    while isfile(outfile)
        outfile = base * "_$(k).jld2"
        k += 1
    end
    return outfile
end

seed = parse_arg(ARGS, 1, 1, Int)
sampling_chains = parse_arg(ARGS, 2, 2, Int)
n_sample = parse_arg(ARGS, 3, 5 * 10^4, Int)
n_layer = parse_arg(ARGS, 4, 10, Int)
show_prog = parse_arg(ARGS, 5, 1, Int) != 0
outdir = length(ARGS) >= 6 ? ARGS[6] : "data/connectance"

c_values = [1.0, 0.8, 0.6, 0.4, 0.2]
S = 8

mkpath(outdir)
outfile = unique_outfile(outdir, seed)
@info "Running connectance sweep" seed outfile

results = Vector{Dict}()
for c in c_values
    @info "Computing connectance c=$c"
    ec = ecosys_config(S; c=c, seed=seed, energetive=true, dissipative=true, positive=false)
    p = generate_model_system(ec)
    Q_range = select_range(p)
    vols_C = volume_range_C(p, Q_range)
    vols_matr = volume_range_EFD(p, :matr, Q_range;
        n_chains=sampling_chains, n_sample=n_sample, n_layer=n_layer,
        show_prog=show_prog, prog_dt=3.0)
    push!(results, Dict(
        :c => c, :S => S, :seed => seed,
        :Qs => Q_range, :vols_C => vols_C,
        :vols_matr => vols_matr, 
        :pm => vols_matr ./ vols_C, 
    ))
end

@save outfile results
@info "Saved" outfile n=length(results)