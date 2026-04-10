""" HPC script: compute feasibility partition volumes across S={2,4,6,8} for a given seed. """

using JLD2

include(joinpath(@__DIR__, "..", "helper.jl"))
include(joinpath(@__DIR__, "..", "partition.jl"))

seed = parse_arg(ARGS, 1, 1, Int)
σsc = parse_arg(ARGS, 2, 1.0, Float64)
d0 = parse_arg(ARGS, 3, 1.0, Float64)
N0 = parse_arg(ARGS, 4, 1.0, Float64)
sampling_chains = parse_arg(ARGS, 5, 2, Int)
n_sample = parse_arg(ARGS, 6, 5 * 10^4, Int)
n_layer = parse_arg(ARGS, 7, 10, Int)
show_prog = parse_arg(ARGS, 8, 1, Int) != 0
outdir = length(ARGS) >= 9 ? ARGS[9] : "data/revision"

mkpath(outdir)
outfile = unique_outfile(outdir, "partition", seed)

@info "Building partition run" seed σsc d0 N0 sampling_chains n_sample n_layer outfile
results = build_seed_partition_run(
    seed;
    σsc=σsc,
    d0=d0,
    N0=N0,
    include_empty=false,
    sampling_chains=sampling_chains,
    n_sample=n_sample,
    n_layer=n_layer,
    save_centers=false,
    show_prog=show_prog,
    prog_dt=3.0,
)
df_system, df_community = collect_partition_tables(results)

@save outfile df_system df_community
@info "Saved partition run" outfile n_system=nrow(df_system) n_community=nrow(df_community)
