"""A full workflow to generate theoretical ecological systems"""

# using .EnerFeas: EnergyConstrProb

# config can be used at different stages, should essentially contain all the possibilities
struct EcosysConfig
    K::Int
    S_type::Symbol
    conne::Float64
    energetive::Bool
    dissipative::Bool
    m_param::Union{Float64, Symbol}
    d_param::Union{Float64, Symbol}
    N0_param::Union{Float64, Symbol}
    ϵ::Float64
    γ::Float64
    seed::Union{Nothing, Int}
end

#export
function ecosys_config(;K::Int, S_type::Symbol, conne::Float64=0.5, energetive::Bool=false, dissipative::Bool=true, m_param=:identical, d_param=:identical, N0_param=:identical, ϵ=0.5, γ=-0.25, seed=nothing)
    (0.0 < conne < 1.0) || error("connectance must be between 0 and 1")
    S_type in [:indiv, :total] || error("S_type must be :indiv or :total")
    m_param in [:identical, :lognormal] || m_param isa Float64 || error("Invalid m_param")
    d_param in [:identical, :lognormal, :allometric] || d_param isa Float64 || error("Invalid d_param")
    N0_param in [:identical, :lognormal, :allometric] || N0_param isa Float64 || error("Invalid N0_param")

    return EcosysConfig(K, S_type, conne, energetive, dissipative, m_param, d_param, N0_param, ϵ, γ, seed)
end


function make_energetive(σ::Matrix{Float64}) #ensuring + > -
    for a in CartesianIndices(σ) 
        i,j = a[1],a[2]
        if σ[i,j] < 0 && σ[j,i] > 0
            x,y = sort(abs.([σ[i,j], σ[j,i]]))
            σ[i,j] = -x; σ[j,i] = +y;
        elseif σ[i,j] > 0 && σ[j,i] < 0
            x,y = sort(abs.([σ[i,j], σ[j,i]]))
            σ[i,j] = +y; σ[j,i] = -x;
        end
    end
    return σ
end

function make_dissipative(σ::Matrix{Float64}, ϵ::Float64=0.5)
    c = minimum(eigvals((σ + σ')/2), init=0.0)
    σ = σ .+ Diagonal(fill((-c + ϵ), size(σ, 1)))
    return σ
end

function impose_connectance(σ::Matrix{Float64}, c::Float64=1.0)
    K = size(σ, 1)
    mask = trues(K, K)
    for i in 1:K, j in 1:K
        if i != j
            mask[i, j] = rand() < c
        end
    end
    return σ .* mask
end

function generate_sigma(K::Int, energetive::Bool, dissipative::Bool, c::Float64=1.0)
    σ = impose_connectance(randn(K,K), c)
    if energetive 
        return make_energetive(σ)
    end
    if dissipative 
        return make_dissipative(σ)
    end
    return σ
end

function generate_bodymass(K::Int, m_param::Union{Float64, Symbol})
    if m_param isa Float64
        return fill(m_param, K)
    elseif m_param == :identical
        return fill(1.0, K)
    elseif m_param == :lognormal
        return exp.(randn(K))
    end
end

# need to change pipeline if further tailering lognormal is needed
function generate_demands(K::Int, d_param::Union{Float64, Symbol}, m::Vector{Float64}, γ::Float64=-0.25)
    if d_param isa Float64
        return fill(d_param, K)  
    elseif d_param == :identical
        return fill(1.0, K)
    elseif d_param == :lognormal
        return exp.(randn(K))
    elseif d_param == :allometric
        return m .^ γ
    end
end

function generate_N0(K::Int, N0_param::Union{Float64, Symbol}, m::Vector{Float64}, γ::Float64=-0.25)
    if N0_param isa Float64
        return fill(N0_param, K)
    elseif N0_param == :identical
        return fill(1.0, K)
    elseif N0_param == :lognormal # we may want to rescale this
        return exp.(randn(K))
    elseif N0_param == :allometric
        return m .^ (1 + γ)
    end
end

#export
function generate_sigma_arrays(ec::EcosysConfig, n_samples::Int=100)
    isnothing(ec.seed) || Random.seed!(ec.seed)
    σs = [generate_sigma(ec.K, ec.energetive, ec.dissipative, ec.conne) for _ in 1:n_samples]
    if n_samples == 1
        σs = σs[1]
    end
    return σs
end

#export
function generate_problem(ec::EcosysConfig, σ::Matrix{Float64})
    m = generate_bodymass(ec.K, ec.m_param)
    d = generate_demands(ec.K, ec.d_param, m, ec.γ)
    N⁰ = generate_N0(ec.K, ec.N0_param, m, ec.γ)
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    S = 0.0
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,ec.K,S,ec.S_type)
end







# function test_proport(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64)
#     size(σ,2) != length(m) && throw(ArgumentError("dimension mismatch between σ and m"))
#     K = size(σ,2); d = fill(0.0, K); N⁰ = fill(0.0, K)
#     Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
#     return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:indiv)
# end

# function test_linear(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64)
#     size(σ,2) != length(m) && throw(ArgumentError("dimension mismatch between σ and m"))
#     K = size(σ,2); d = fill(1.0, K); N⁰ = fill(1.0, K); 
#     Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
#     return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:indiv)
# end

# function test_quadratic(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64, ϵ::Float64=0.8)
#     size(σ,2) != length(m) && throw(ArgumentError("dimension mismatch between σ and m"))
#     K = size(σ, 2)
#     c = minimum(eigvals((σ + σ')/2), init=0.0)
#     σ += Diagonal(fill((-c + ϵ), K))
#     d = fill(0.1, K); N⁰ = fill(0.2, K); # different critical values for d and N⁰
#     Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
#     return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:total)
# end

# function test_linear_compare(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64, ϵ::Float64=0.8)
#     size(σ,2) != length(m) && throw(ArgumentError("dimension mismatch between σ and m"))
#     K = size(σ, 2)
#     c = minimum(eigvals((σ + σ')/2), init=0.0)
#     σ += Diagonal(fill((-c + ϵ), K))
#     d = fill(0.1, K); N⁰ = fill(0.2, K); # different critical values for d and N⁰
#     Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
#     return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:linear)
# end

# function test_allom(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64, γ::Float64=-0.25)
#     size(σ,2) != length(m) && throw(ArgumentError("dimension mismatch between σ and m"))
#     K = size(σ, 2); d = m .^ γ; N⁰ = m .^ (1+γ)
#     Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
#     return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:indiv)
# end


