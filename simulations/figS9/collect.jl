""" Collect connectance sweep results across seeds and export to JSON. """

using JLD2
using JSON

datadir = length(ARGS) >= 1 ? ARGS[1] : "data/connectance"
outfile = length(ARGS) >= 2 ? ARGS[2] : "data/connectance/collected.json"

rows = []
for f in readdir(datadir; join=true)
    endswith(f, ".jld2") || continue
    results = load(f, "results")
    for d in results
        push!(rows, Dict(
            "seed" => d[:seed],
            "c"    => d[:c],
            "Qs"   => d[:Qs],
            "pm"   => d[:pm],
        ))
    end
end

@info "Collected $(length(rows)) rows from $(datadir)"
open(outfile, "w") do io
    JSON.print(io, rows)
end
@info "Saved to $(outfile)"
