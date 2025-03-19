include("src/EnerFeas.jl");
using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots
using DataFrames
using ProgressMeter

include("src/trophic.jl");

function test_proport(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64)
    d = fill(0.0, 3); N⁰ = fill(0.0, 3); K = 3
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:Individual)
end

function test_linear(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64)
    d = fill(0.0, 3); N⁰ = fill(0.2, 3); K = 3
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:Individual)
end

function test_quadratic(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64, ϵ::Float64=0.8)
    c = minimum(eigvals((σ + σ')/2), init=0.0)
    σ += Diagonal(fill((-c + ϵ), 3))
    d = fill(0.0, 3); N⁰ = fill(0.25, 3); K = 3
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:Total)
end

#1 --------- Test proportional domains --------- 
df_plain = DataFrame(rep_id=Int[], type=String[], S_ind=Float64[], norm_vol=Float64[])

S_levels = 0.5:0.025:10.0;
σs = [randn(3, 3) for _ in 1:50]
m0 = fill(1.0, 3)
for i in eachindex(σs)
    σ = σs[i]
    for S_ind in S_levels
        p = test_proport(σ, m0, S_ind)
        vol_d, vol_sim = volume_EFD(p, true)
        norm_vol = vol_d / vol_sim
        push!(df_plain, (i, "proport", S_ind, norm_vol))
        @info "rep_id:$i, S_ind:$S_ind"
    end
end

# plot: for each rep_id, plot the norm_vol - S_ind curve
plot(df_plain.S_ind, df_plain.norm_vol, group=df_plain.rep_id,
    xlabel="S_ind", ylabel="norm_vol", title="Proportional",
    legend=:topright, label=false, markersize=2, color = :grey)



#2 --------- Comparision inside trophic chains ---------
Random.seed!(1010)
tr_types = ["hierarchy", "omnivory", "exploitative", "competitive"]
σs = Dict(tr_type => [generate_trophic(tr_type) for _ in 1:50] for tr_type in tr_types)

#2.1 --------- Test linear domains --------- 
df_linear = DataFrame(rep_id=Int[], tr_type=String[], S=Float64[], norm_vol=Float64[], vol_d=Float64[], vol_sim=Float64[])
S_levels = 0.0:0.05:5.0;

for tr_type in tr_types
    for i in range(1,50)
        σ = σs[tr_type][i]
        pbar = Progress(length(S_levels); desc="tr_type:$tr_type, rep_id:$i")
        for S_ind in S_levels
            p = test_linear(σ, fill(1.0, 3), S_ind)
            vol_d, vol_sim = volume_EFD(p, true) # use exact solu, fast in low dims
            norm_vol = vol_d / vol_sim
            push!(df_linear, (i, tr_type, S_ind, norm_vol, vol_d, vol_sim))
            next!(pbar)
        end
    end
end

#2.2 --------- Test quadratic domains --------- 
df_quadratic = DataFrame(rep_id=Int[], tr_type=String[], S=Float64[], norm_vol = Float64[], vol_d=Float64[], vol_sp=Float64[])
S_levels = 0.0:0.1:5.0;

for tr_type in tr_types
    for i in range(1,50) # only use half
        σ = σs[tr_type][i]
        pbar = Progress(length(S_levels); desc="tr_type:$tr_type, rep_id:$i")
        for S_tot in S_levels
            p = test_quadratic(σ, fill(1.0, 3), S_tot)
            # @info "tr_type:$tr_type, rep_id:$i, S_tot:$S_tot"
            vol_d, vol_sp = volume_EFD(p, false, 5,1)
            norm_vol = vol_d / vol_sp
            push!(df_quadratic, (i, tr_type, S_tot, norm_vol, vol_d, vol_sp))
            next!(pbar)
        end
        # @info "tr_type:$tr_type, rep_id:$i"
    end
end

using CSV
CSV.write("analysis/trophic_data_linear.csv", df_linear)
CSV.write("analysis/trophic_data_quadratic.csv", df_quadratic)