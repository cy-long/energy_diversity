include("src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots, StatsPlots
using DataFrames
using ProgressMeter
using IterTools



#2 --------- Comparision inside trophic chains ---------
Random.seed!(1010)
tr_types = ["chain", "omnivory", "exploitative", "competitive"]
σs = Dict(
    tr_type => [generate_trophic(tr_type) for _ in 1:50] for tr_type in tr_types
)
cm = Dict(
    "chain" => :red, "omnivory" => :blue, "exploitative" => :green, "competitive" => :purple
)

#2.1 --------- Test linear domains --------- 
df_linear = DataFrame(rep_id=Int[], tr_type=String[], S=Float64[], norm_vol=Float64[], vol_d=Float64[], vol_sim=Float64[])
S_levels = 0.0:0.05:20.0;

for (i, tr_type) in Iterators.product(1:50, tr_types)
    σ = σs[tr_type][i]
    m = sort(rand(LogNormal(0.25, 1.0),3),rev=true)
    # m = fill(1.0, 3) # should recover previous
    for S_ind in S_levels
        # p = test_linear(σ, fill(1.0, 3), S_ind)
        p = test_allom(σ, m,  S_ind, -0.25) # use allometric scaling
        vol_d, vol_sim = volume_EFD(p, true) # use exact solu, fast in low dims
        norm_vol = vol_d / S^3
        push!(df_linear, (i, tr_type, S_ind, norm_vol, vol_d, vol_sim))
    end
    @info "tr_type:$tr_type, rep_id:$i"
end


@df df_linear plot(:S, :norm_vol, 
    group=:tr_type,
    facet=:tr_type,
    color=[cm[t] for t in df_linear.tr_type],
    layout=(2,2),
    xlabel="S_ind", ylabel="norm_vol",
    ylim=(0,1.0),
    markersize=2)

plot_path = "figures/linear_.pdf"
plot!(size=(800, 600), titlefont_size=16, legendfont_size=12)
savefig(plot_path)

using CSV
CSV.write("analysis/trophic_linear_am.csv", df_linear)


# ps_linear = DataFrame(rep_id=Int[], tr_type=String[], Sc=Float64[], Ω=Float64[])
# for (type, i) in Iterators.product(tr_types, 1:50)
#     dfi = df_linear[(df_linear.tr_type .== type) .& (df_linear.rep_id .== i), :]
#     dfi.norm_vol[isnan.(dfi.norm_vol)] .= 0.0 # treat NaN as 0
#     Sc = maximum(dfi.S[dfi.vol_d .== 0.0])
#     Ω = maximum(dfi.norm_vol)
#     push!(ps_linear, (i, type, Sc, Ω))
# end
# using CSV
# CSV.write("analysis/trophic_data_linear.csv", df_linear)
# CSV.write("analysis/omega_sc_linear.csv", ps_linear)

# scatter(ps_linear.Sc, ps_linear.Ω, 
#     group=ps_linear.tr_type,
#     xlabel="S_c", ylabel="Ω",
#     legend=:topright, markersize=5,
#     xlim = (0,5),
#     color=[cm[t] for t in ps_linear.tr_type])


# #2.2 --------- Test quadratic domains --------- 


# df_quadratic = DataFrame(rep_id=Int[], tr_type=String[], S=Float64[], norm_vol = Float64[], vol_d=Float64[], vol_sp=Float64[])
# S_levels = 0.0:0.1:5.0;

# for tr_type in tr_types
#     for i in range(1,50) # only use half
#         σ = σs[tr_type][i]
#         pbar = Progress(length(S_levels); desc="tr_type:$tr_type, rep_id:$i")
#         for S_tot in S_levels
#             p = test_quadratic(σ, fill(1.0, 3), S_tot)
#             # @info "tr_type:$tr_type, rep_id:$i, S_tot:$S_tot"
#             vol_d, vol_sp = volume_EFD(p, false, 5,1)
#             norm_vol = vol_d / vol_sp
#             push!(df_quadratic, (i, tr_type, S_tot, norm_vol, vol_d, vol_sp))
#             next!(pbar)
#         end
#         # @info "tr_type:$tr_type, rep_id:$i"
#     end
# end

# # df_quadratic.vol_d = replace(df_quadratic.vol_d, NaN => 0.0)
# color_values = [cm[t] for t in df_quadratic.tr_type]
# @df df_quadratic plot(:S, df_quadratic.vol_d ./ df_quadratic.S.^3, 
#         group=:tr_type,
#         facet=:tr_type,
#         color=color_values, # Explicitly match color with cleaned data
#         layout=(2,2),
#         xlabel="S_tot", ylabel="rescaled vol.",
#         ylim=(0,0.3),
#         markersize=2,
#         guide_position=:outer,
#         plot_title="Quadratic domains; d=0.1, N⁰=0.2")

# # save plot into pdf
# plot_path = "figures/quad_1.pdf"
# plot!(size=(800, 600), titlefont_size=16, legendfont_size=12)
# savefig(plot_path)



# using CSV
# CSV.write("analysis/trophic_data_linear.csv", df_linear)
# CSV.write("analysis/trophic_data_quadratic_dN.csv", df_quadratic)


# #1 --------- Test proportional domains --------- 


# df_plain = DataFrame(rep_id=Int[], type=String[], S_ind=Float64[], norm_vol=Float64[])

# S_levels = 0.5:0.025:10.0;
# σs = [randn(3, 3) for _ in 1:50]
# m0 = fill(1.0, 3)
# for i in eachindex(σs)
#     σ = σs[i]
#     for S_ind in S_levels
#         p = test_proport(σ, m0, S_ind)
#         vol_d, vol_sim = volume_EFD(p, true)
#         norm_vol = vol_d / vol_sim
#         push!(df_plain, (i, "proport", S_ind, norm_vol))
#         @info "rep_id:$i, S_ind:$S_ind"
#     end
# end

# plot: for each rep_id, plot the norm_vol - S_ind curve
# plot(df_plain.S_ind, df_plain.norm_vol, group=df_plain.rep_id,
#     xlabel="S_ind", ylabel="norm_vol", title="Proportional",
#     legend=:topright, label=false, markersize=2, color = :grey)



#4 ---- casecade estimation ---- Will that amplify the error?
