""" Testing fast volume calculation for series of quadratic domains """

include("../src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots, StatsPlots
using DataFrames, ProgressMeter, IterTools

S_range = Vector(1e-6:0.1:5.0);

# ec0 = ecosys_config(K=4, S_type=:total, conne = 0.5, m_param=:lognormal, d_param=:allometric, N0_param=:allometric, seed=25);
ec0 = ecosys_config(K=4, S_type=:total, conne = 0.5, d_param=0.08, N0_param=0.1, seed=100);
σs = generate_sigma_arrays(ec0, 25);

df0 = DataFrame(id=Int[], S=Float64[], res_vol=Float64[]);

for i in eachindex(σs)
    @info "rep_id:$i"
    p0 = generate_problem(ec0, σs[i]);
    vol0 = volume_cascade_EFD(p0, S_range);
    for (j, S) in enumerate(S_range)
        push!(df0, (i, S, vol0[j] / S^ec0.K))
    end
end

@df df0 plot(:S, :res_vol, 
    group=:id,
    color=:darkorange,
    label=nothing,
    xlabel="S_tot", ylabel="res. vol",
    ylim = (0,2.5),
    markersize=2)