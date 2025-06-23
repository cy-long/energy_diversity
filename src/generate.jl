"""The toolset to generate theoretical ecological systems"""

struct EcosysConfig
    K::Int
    S_type::Symbol
    n_scale::Float64
    n_var::Float64
    conne::Float64
    energetive::Bool
    dissipative::Bool
    k_param::Union{Float64, Symbol}
    d_param::Union{Float64, Symbol}
    N0_param::Union{Float64, Symbol}
    ϵ::Float64
    γ::Float64
    seed::Union{Nothing, Int}
end

function ecosys_config(;K::Int, S_type::Symbol, n_scale::Float64=1.0, n_var::Float64=1.0, conne::Float64=1.0, energetive::Bool=false, dissipative::Bool=true, k_param=:identical, d_param=:identical, N0_param=:identical, ϵ=0.5, γ=-0.25, seed=nothing)
    (0.0 < conne <= 1.0) || error("connectance must be between 0 and 1")
    S_type in [:indiv, :total] || error("S_type must be :indiv or :total")
    k_param in [:unit, :identical, :normal] || k_param isa Float64 || error("Invalid k_param")
    d_param in [:identical, :lognormal, :allometric] || d_param isa Float64 || error("Invalid d_param")
    N0_param in [:identical, :lognormal] || N0_param isa Float64 || error("Invalid N0_param")

    return EcosysConfig(K, S_type, n_scale, n_var, conne, energetive, dissipative, k_param, d_param, N0_param, ϵ, γ, seed)
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

function make_dissipative!(σ::Matrix{Float64}, ϵ::Float64=0.25)
    c = minimum(eigvals((σ + σ') / 2), init=0.0)
    shift = -c + ϵ
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

function generate_sigma(K::Int, energetive::Bool, dissipative::Bool, n_scale::Float64=1.0, n_var::Float64=1.0, c::Float64=1.0, ϵ::Float64=0.25)
    σ = rand(Normal(0, n_var),K,K)
    if c < 1.0
        impose_connectance!(σ, c)
    end
    if energetive
        make_energetive!(σ)
    end
    if dissipative 
        make_dissipative!(σ, ϵ)
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

function generate_k(K::Int, k_param::Union{Float64, Symbol})
    if k_param isa Float64
        return fill(k_param, K)
    elseif k_param == :lognormal
        return exp.(randn(K))
    elseif k_param == :unit
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
    σs = [generate_sigma(ec.K, ec.energetive, ec.dissipative, ec.n_scale, ec.n_var, ec.conne) for _ in 1:n_samples]
    if n_samples == 1
        σs = σs[1]
    end
    return σs
end

function generate_problem(ec::EcosysConfig, σ::Matrix{Float64})
    k = generate_k(ec.K, ec.k_param)
    N⁰ = generate_N0(ec.K, ec.N0_param, fill(0.0, ec.K), ec.γ)
    d = generate_d(ec.K, ec.d_param, N⁰, ec.γ)
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    S = 0.0
    return EnergyConstrProb(σ,Λ,Q,c,k,d,N⁰,ec.K,S,ec.S_type)
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