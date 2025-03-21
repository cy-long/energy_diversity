using Revise

includet("src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots, StatsPlots
using DataFrames
using ProgressMeter
using IterTools


Random.seed!(1101)
σ =  generate_trophic("exploitative");
m = rand(LogNormal(0.25, 0.1), 3);
p0 = test_quadratic(σ, m, 5.0);
S_range = Vector(0.1:0.05:5.0);

vol1 = volume_cascade_EFD(p0, S_range, 1.0);
vol2 = volume_cascade_EFD(p0, S_range, 0.5);


# @elapsed vol2 = [first(volume_EFD(test_quadratic(σ, m, S))) for S in S_range]

plot(xlabel="S", ylabel="volume", title="Volume of EFD vs S", legend=:topright, markersize=2, grid=false);
plot!(S_range, vol1./ S_range .^3, label="cascade", linewidth=1)
plot!(S_range, vol2./ S_range .^3, label="cascade+ball", linewidth=2)