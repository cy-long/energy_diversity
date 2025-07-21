"""The toolset to generate theoretical ecological systems"""

struct EcosysConfig
    K::Int
    S_type::Symbol
    n_scale::Float64
    n_var::Float64
    conne::Float64
    energetive::Bool
    dissipative::Bool
    positive::Bool
    ϵ_param::Union{Float64, Symbol}
    d_param::Union{Float64, Symbol}
    N0_param::Union{Float64, Symbol}
    σ0::Float64
    γ::Float64
    seed::Union{Nothing, Int}
end

function generate_model_system(size::Int, type::Symbol, seed::Int, σ_scale::Float64, d0::Float64, N0::Float64)
    Random.seed!(seed)
    type ∈ [:total, :indiv] || error("Type must be either :total or :indiv")
    σ = generate_sigma(size, true, true, false, σ_scale, 1.0, 1.0, 0.5);
    d = d0 * rand(Normal(1.0, 0.1), size)
    N⁰ = N0 * rand(Normal(1.0, 0.1), size)
    Λ = inv(σ); Q = (Λ + Λ') / 2; c = -Λ * d;
    ϵ = fill(0.0, size);
    p = EnergyConstrProb(σ, Λ, Q, c, ϵ, d, N⁰, size, 0.0, type);
    return p
end

function sub_model_system(size::Int, p::EnergyConstrProb)
    size ≤ p.K || error("Size must be less than or equal to the number of populations in the model system")
    if size == p.K
        return p
    end
    σ = p.σ[1:size, 1:size]
    Λ = inv(σ); Q = (Λ + Λ') / 2; c = -Λ * p.d[1:size];
    return EnergyConstrProb(σ, Λ, Q, c, fill(0.0, size), p.d[1:size], p.N⁰[1:size], size, 0.0, p.type)
end

function volume_range_flux(p::EnergyConstrProb, Q_range::Vector{Float64})
    return [Q^p.K/(prod(p.N⁰)*factorial(p.K)) for Q in Q_range]
end

# function select_range(p::EnergyConstrProb)
#     Q1 = baseline_supply(p);
#     Q2 = optimal_supply(p);
#     Q_range = vcat(
#         range(1e-4, Q1; length=25)[1:end], # scan the critical
#         range(Q1, 2*Q2; length=350+1)[2:end], # scan the optimum
#         range(2*Q2, max(10*Q2, 100.0); length=125+1)[2:end], # at least 10^2
#     )
#     return Q_range
# end

function select_range(p::EnergyConstrProb)
    Q1 = baseline_supply(p)
    Q2 = optimal_supply(p)
    Q_range = vcat(
        exp.(range(log(1e-4), log(Q1); length=25)),            # before critical
        exp.(range(log(Q1), log(Q2); length=351)[2:end]),     # around optimum
        exp.(range(log(Q2), log(max(10Q2, 100.0)); length=126)[2:end])  # tail
    )
    return Q_range
end

# ---- old functions, need to be updated ----

# Reconstruct this function accroding to new texts
function ecosys_config(;
    K::Int, S_type::Symbol, n_scale::Float64=1.0, n_var::Float64=1.0, conne::Float64=1.0, 
    energetive::Bool=false, dissipative::Bool=true, positive::Bool=false,
    d_param=:identical, d_scale=1.0,
    N0_param=:identical, N0_scale=1.0,
    ϵ_param=0.0, ϵ_scale=0.0,
    σ0=0.25, γ=-0.25, seed=nothing
)

    (0.0 < conne <= 1.0) || error("connectance must be between 0 and 1")
    S_type in [:indiv, :total] || error("S_type must be :indiv or :total")
    ϵ_param in [:unit, :identical, :normal] || ϵ_param isa Float64 || error("Invalid ϵ_param")
    d_param in [:identical, :lognormal, :allometric] || d_param isa Float64 || error("Invalid d_param")
    N0_param in [:identical, :lognormal] || N0_param isa Float64 || error("Invalid N0_param")

    return EcosysConfig(K, S_type, n_scale, n_var, conne, energetive, dissipative, positive, ϵ_param, d_param, N0_param, σ0, γ, seed)
end

function make_energetive!(σ::Matrix{Float64})  # modifies σ in-place
    for a in CartesianIndices(σ)
        i, j = a[1], a[2]
        if σ[i,j] < 0 && σ[j,i] > 0
            x, y = sort(abs.([σ[i,j], σ[j,i]]))
            σ[i,j] = -x
            σ[j,i] =  y
        elseif σ[i,j] > 0 && σ[j,i] < 0
            x, y = sort(abs.([σ[i,j], σ[j,i]]))
            σ[i,j] =  y
            σ[j,i] = -x
        end
    end
    return nothing
end

function make_dissipative!(σ::Matrix{Float64}, σ0::Float64=0.25)
    c = minimum(eigvals((σ + σ') / 2), init=0.0)
    shift = -c + σ0
    for i in 1:size(σ, 1)
        σ[i, i] += shift
    end
    return nothing
end

function impose_connectance!(σ::Matrix{Float64}, c::Float64=1.0)
    K = size(σ, 1)
    for i in 1:K, j in 1:K
        if i != j && rand() ≥ c
            σ[i,j] = 0.0
        end
    end
    return nothing
end

function generate_sigma(K::Int, energetive::Bool, dissipative::Bool, positive::Bool, n_scale::Float64=1.0, n_var::Float64=1.0, c::Float64=1.0, σ0::Float64=0.25)
    σ = rand(Normal(0, n_var),K,K)
    if c < 1.0
        impose_connectance!(σ, c)
    end
    if positive
        σ = abs.(σ)
    end
    if energetive
        make_energetive!(σ)
    end
    if dissipative 
        make_dissipative!(σ, σ0)
    end
    return n_scale * σ
end

# need to change pipeline if lognormal is needed
function generate_d(K::Int, d_param::Union{Float64, Symbol}, m::Vector{Float64}, γ::Float64=-0.25)
    if d_param isa Float64
        return fill(d_param, K)  
    elseif d_param == :identical
        return fill(1.0, K)
    elseif d_param == :lognormal
        return 0.2*exp.(randn(K))
    elseif d_param == :allometric
        return m .^ γ
    end
end

function generate_eps(K::Int, ϵ_param::Union{Float64, Symbol})
    if ϵ_param isa Float64
        return fill(ϵ_param, K)
    elseif ϵ_param == :lognormal
        return exp.(randn(K))
    elseif ϵ_param == :unit
        return fill(1.0, K)
    end
end

function generate_N0(K::Int, N0_param::Union{Float64, Symbol}, m::Vector{Float64}, γ::Float64=-0.25)
    if N0_param isa Float64
        return fill(N0_param, K)
    elseif N0_param == :identical
        return fill(1.0, K)
    elseif N0_param == :lognormal # we may want to rescale this
        return rand(LogNormal(-0.5, 1.0), K)
    end
end

function generate_sigma_arrays(ec::EcosysConfig, n_samples::Int=100)
    isnothing(ec.seed) || Random.seed!(ec.seed)
    σs = [generate_sigma(ec.K, ec.energetive, ec.dissipative, ec.positive, ec.n_scale, ec.n_var, ec.conne, ec.σ0) for _ in 1:n_samples]
    if n_samples == 1
        σs = σs[1]
    end
    return σs
end

function generate_problem(ec::EcosysConfig, σ::Matrix{Float64})
    ϵ = generate_eps(ec.K, ec.ϵ_param)
    N⁰ = generate_N0(ec.K, ec.N0_param, fill(0.0, ec.K), ec.γ)
    d = generate_d(ec.K, ec.d_param, N⁰, ec.γ)
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    S = 0.0
    return EnergyConstrProb(σ,Λ,Q,c,ϵ,d,N⁰,ec.K,S,ec.S_type)
end


# Tried to analyze trophic systems but didn't find immediate results.
# function generate_trophic_chain(K::Int, n_scale::Float64=1.0, n_var::Float64=1.0, effi_mean::Float64=0.5, effi_var::Float64=0.1)
#     σ = abs.(Matrix(Diagonal(rand(Normal(0, n_var),K))))
#     uptake = n_scale * abs.(rand(Normal(0, n_var), K-1))
#     effi = rand(Normal(effi_mean, effi_var), K-1); effi[effi .< 0.0] .= 0.0
#     intake = -effi .* uptake
#     for i in 1:K-1
#         σ[i, i+1] = uptake[i]
#         σ[i+1, i] = intake[i]
#     end
#     make_dissipative!(σ)
#     return σ
# end


# function generate_simple_chain(K::Int, a::Float64=2.0, b::Float64=1.0, c::Float64=1.0)
#     σ = diagm(a * ones(K))
#     for i in 1:K-1
#         σ[i, i+1] = c
#         σ[i+1, i] = -b
#     end
#     make_dissipative!(σ)
#     return σ
# end


# f: fraction of +-, --, ++; δ: probability of perturbation for +-;
# signs are reversed in σ matrix; connectance set as 1.0;
# setup a cascade-like trophic chain, lower # eaten by higher #
# function generate_trophic(K::Int, f_np::Float64=1.0, f_nn::Float64=0.0, f_pp::Float64=0.0)
#     if f_np + f_nn + f_pp > 1.0 || f_np < 0.0 || f_nn < 0.0 || f_pp < 0.0
#         throw(ArgumentError(("invalid distribution of sign patterns")))
#     end
#     σ = abs.(randn(K, K))
#     signs = zeros(Int8, K, K)
#     for i in 1:K-1, j in i+1:K
#         r = rand()
#         if r < f_pp
#             signs[i,j] = 1; signs[j,i] = 1
#         elseif r < f_pp + f_nn
#             signs[i,j] = -1; signs[j,i] = -1
#         elseif r < f_pp + f_nn + f_np
#             signs[i,j] = -1; signs[j,i] = 1
#             x,y = sort(abs.([σ[i,j], σ[j,i]]))
#             σ[i,j] = y; σ[j,i] =x
#         end
#     end
#     σ = -σ .* signs # reverse the sign because of the convention
#     make_dissipative!(σ)
#     return σ
# end


# function perturb_trophic(σ::AbstractMatrix, δ::Float64=0.0)
#     σ_new = copy(σ)
#     K = size(σ, 1)
#     for i in 1:K-1, j in i+1:K
#         if σ_new[i,j]*σ_new[j,i] < 0.0 && rand() < δ # variance might be high
#             σ_new[i,j] = -σ_new[i,j]
#             σ_new[j,i] = -σ_new[j,i]
#         end
#     end
#     make_dissipative!(σ_new, 0.1)
#     return σ_new
# end