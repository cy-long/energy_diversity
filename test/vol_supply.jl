""" Fast volume calculation for a range of total total_supply (ploting the curves) """

include("../src/EnerFeas.jl");

using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots, StatsPlots
using DataFrames, ProgressMeter, IterTools

S_range = Vector(1e-6:0.1:5.0);

ec = ecosys_config(K=4, S_type=:total, conne = 0.5, d_param=0.08, N0_param=0.1, seed=100);
σs = generate_sigma_arrays(ec, 5);

df = DataFrame(id=Int[], S=Float64[], res_vol=Float64[]);
for i in eachindex(σs)
    @info "rep_id:$i"
    p = generate_problem(ec, σs[i]);
    vol = volume_cascade_EFD(p, S_range);
    for (j, S) in enumerate(S_range)
        push!(df, (i, S, vol[j] / S^ec.K))
    end
end

@df df plot(:S, :res_vol, 
    group=:id,
    # color=:dodgerblue,
    label=nothing,
    xlabel="S_tot", ylabel="res. vol",
    ylim = (0,1.0),
    markersize=2)

# Linear counterparts
ec1 = ecosys_config(K=4, S_type=:indiv, conne = 0.5, d_param=0.08, N0_param=0.1, seed=100);
σs1 = generate_sigma_arrays(ec1, 5); # same sigma as before

df1 = DataFrame(id=Int[], S=Float64[], res_vol=Float64[]);
for i in eachindex(σs)
    @info "rep_id:$i"
    p = generate_problem(ec1, σs1[i]);
    for (j, S) in enumerate(S_range)
        p.S = S
        vol_j, vol_all = volume_EFD(p, true);
        push!(df1, (i, S, vol_j/((ec1.K*S)^ec1.K)))
    end
end
df1[isnan.(df1.res_vol), :res_vol] .= 0.0;


@df df1 plot(:S, :res_vol, 
    group=:id,
    label=nothing,
    xlabel="S_ind", ylabel="res. vol",
    ylim = (0,0.1),
    markersize=2)