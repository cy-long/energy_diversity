"""Encode specific network"""

# rms = [randn(10,10) for _ in 1:100000]
# eigs = [eigvals((a+a')/2) for a in rms]
# eigmins = [minimum(e) for e in eigs]
# histogram(eigmins, bins=100, xlabel="Eigenvalues", ylabel="Frequency", title="Min Eigvals of random symmetric mats")

using .EnerFeas: EnergyConstrProb

function purely_hierarchy()
    σ = abs.(randn(3,3)) 
    if σ[1,2] < σ[2,1] 
        σ[2,1], σ[1,2] = sort!(abs.(rand(2)))
    end
    if σ[2,3] < σ[3,2] 
        σ[3,2], σ[2,3] = sort!(abs.(rand(2)))
    end
    return -σ .*[-1 -1 0; +1 -1 -1; 0 +1 -1]
end

function omnivory()
    σ = abs.(randn(3,3))
    if σ[1,2] < σ[2,1] 
        σ[2,1], σ[1,2] = sort!(abs.(rand(2)))
    end
    if σ[2,3] < σ[3,2] 
        σ[3,2], σ[2,3] = sort!(abs.(rand(2)))
    end
    if σ[1,3] < σ[3,1]
        σ[3,1], σ[1,3] = sort!(abs.(rand(2)))
    end
    return -σ .*[-1 -1 -1; +1 -1 -1; +1 +1 -1]
end

# competition within the same trophic level are temporarily omitted
function exploitative()
    σ = abs.(randn(3,3))
    if σ[1,2] < σ[2,1] 
        σ[2,1], σ[1,2] = sort!(abs.(rand(2)))
    end
    if σ[1,3] < σ[3,1]
        σ[3,1], σ[1,3] = sort!(abs.(rand(2)))
    end
    return -σ .*[-1 -1 -1; +1 -1 0; +1 0 -1]
end

function purely_competitive()
    σ = abs.(randn(3,3))
    return σ
end

function generate_trophic(type::String)
    if type == "chain"
        σ = purely_hierarchy()
    elseif type == "omnivory"
        σ = omnivory()
    elseif type == "exploitative"
        σ = exploitative()
    elseif type == "competitive"
        σ = purely_competitive()
    else
        throw(ArgumentError("trophic types not accepted"))
    end
    return σ
end


# These "test_" functions are truly an interface from σ, d, m to the EnergyConstrProb!
# They can only be used for K=3 now
# Use them utill we update the interface properly (plan to give interface for adding customized constraints)
function test_proport(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64)
    K = 3; d = fill(0.0, K); N⁰ = fill(0.0, K)
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:indiv)
end

function test_linear(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64)
    K = 3; d = fill(1.0, K); N⁰ = fill(1.0, K); 
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:indiv)
end

function test_quadratic(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64, ϵ::Float64=0.8)
    c = minimum(eigvals((σ + σ')/2), init=0.0)
    σ += Diagonal(fill((-c + ϵ), 3))
    
    K = 3; d = fill(0.1, K); N⁰ = fill(0.2, K); # different critical values for d and N⁰
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:total)
end

function test_allom(σ::Matrix{Float64}, m::Vector{Float64}, S::Float64, γ::Float64=-0.25)
    K = 3; d = m .^ γ; N⁰ = m .^ (1+γ)
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d
    return EnergyConstrProb(σ,Λ,Q,c,m,d,N⁰,K,S,:indiv)
end


# Calculate the baseline supply needed to sustain the ecosystem with least biomass (averaged by K)
baseline_supply(p::EnergyConstrProb) = sum(p.d + inv(p.Λ) * p.N⁰) / p.K

# Calculate the individual supply at state s (averaged by K)
individual_supply(s::Vector{Float64}, p::EnergyConstrProb) = dot(s, p.m) / p.K

# Calculate the total supply at state s (weighed by N)
total_supply(s::Vector{Float64}, p::EnergyConstrProb) = transpose(s) * p.Λ * (s - p.d)