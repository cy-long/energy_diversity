"""The toolset to generate model ecological systems"""
struct EcosysConfig
    S::Int
    σ_scale::Float64
    σ0::Float64
    energetive::Bool
    dissipative::Bool
    positive::Bool
    ϵ::Union{Float64, Symbol}
    d0::Union{Float64, Symbol}
    N0::Union{Float64, Symbol}
    seed::Union{Nothing, Int}
end

function ecosys_config(S::Int;
    σ_scale::Float64=1.0, σ0=0.25,
    energetive::Bool=false, dissipative::Bool=true, positive::Bool=false,
    d0=1.0, N0=1.0, ϵ=0.0, seed=nothing
)
    ϵ ≥ 0.0 || error("Invalid ϵ_param")
    d0 > 0.0 || error("Invalid d_param")
    N0 > 0.0 || error("Invalid N0_param")
    σ_scale > 0.0 || error("Invalid σ_scale")
    σ0 ≥ 0.0 || error("Invalid σ0")
    return EcosysConfig(S, σ_scale, σ0, energetive, dissipative, positive, ϵ, d0, N0, seed)
end

function generate_sigma_arrays(ec::EcosysConfig, n_samples::Int=100)
    isnothing(ec.seed) || Random.seed!(ec.seed)
    σs = [generate_sigma(ec.S, ec.energetive, ec.dissipative, ec.positive, ec.σ_scale, 1.0, ec.σ0) for _ in 1:n_samples]
    if n_samples == 1
        σs = σs[1]
    end
    return σs
end

function make_energetive!(σ::Matrix{Float64})
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

function make_dissipative!(σ::Matrix{Float64}, σ0::Float64)
    c = minimum(eigvals((σ + σ') / 2), init=0.0)
    shift = -c + σ0
    for i in 1:size(σ, 1)
        σ[i, i] += shift
    end
    return nothing
end

function impose_connectance!(σ::Matrix{Float64}, c::Float64)
    K = S(σ, 1)
    for i in 1:K, j in 1:K
        if i != j && rand() ≥ c
            σ[i,j] = 0.0
        end
    end
    return nothing
end

function generate_sigma(
    S::Int, energetive::Bool, dissipative::Bool, positive::Bool,
    σ_scale::Float64=1.0, c::Float64=1.0, σ0::Float64=0.5
    )
    σ = rand(Normal(0, 1.0),S,S)
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
    return σ_scale * σ
end

function generate_sigma_arrays(S, seed, σ_sc, σ0, n::Int)::Vector{Matrix{Float64}}
    σs = Vector{Matrix{Float64}}(undef, n)
    Random.seed!(seed)
    for i in 1:n
        σs[i] = generate_sigma(
            S, false, true, false,
            σ_sc, 1.0, σ0
        )
    end
    return σs
end

function generate_model_system(
    S::Int, seed::Int, σ_scale::Float64, d0::Float64, N0::Float64;
    d_sd::Float64=0.1, N0_sd::Float64=0.1, σ0::Float64=0.5
)::EcosysParams
    Random.seed!(seed)
    σ = generate_sigma(S, true, true, false, σ_scale, 1.0, σ0);
    d = d0 * rand(Normal(1.0, d_sd), S)
    N⁰ = N0 * rand(Normal(1.0, N0_sd), S)
    ϵ = zeros(S);
    p = EcosysParams(S, σ, d, 0.0, N⁰, ϵ);
    return p
end

function generate_model_system(ec::EcosysConfig, σ::Matrix{Float64})::EcosysParams
    ϵ = zeros(ec.S)
    N⁰ = fill(ec.N0, ec.S)
    d = fill(ec.d0, ec.S)
    return EcosysParams(ec.S, σ, d, 0.0, N⁰, ϵ)
end

function sub_model_system(S::Int, p::EcosysParams)::EcosysParams
    S ≤ p.S || error("Size exceeds the original model")
    if S == p.S
        return p
    end
    return EcosysParams(S, p.σ[1:S, 1:S], p.d[1:S], 0.0, p.N⁰[1:S], zeros(S))
end