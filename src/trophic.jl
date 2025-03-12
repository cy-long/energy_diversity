"""Encode specific network"""

using Random
using Distributions
using LinearAlgebra
include("_problem.jl")

function purely_hierarchy()
    σ = abs.(randn(3,3)) 
    if σ[1,2] < σ[2,1] 
        σ[2,1], σ[1,2] = sort!(abs.(rand(2)))
    end
    if σ[2,3] < σ[3,2] 
        σ[3,2], σ[2,3] = sort!(abs.(rand(2)))
    end
    return(-σ .*[-1 -1 0; +1 -1 -1; 0 +1 -1])
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
    return(-σ .*[-1 -1 -1; +1 -1 -1; +1 +1 -1])
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
    return(-σ .*[-1 -1 -1; +1 -1 0; +1 0 -1])
end


function purely_competitive()
    σ = abs.(randn(3,3))
    return(σ)
end

function generate_trophic(type::String)
    if type == "hierarchy"
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
    return(σ)
end


# Let's consider the simplest case: same m,d,N⁰
function create_ep_trophic(σ::Matrix{Float64}, S::Float64, type::Symbol)
    if type ∉ [:Individual, :Total]
        throw(ArgumentError("type must be either :Individual or :Total"))
    end
    m = fill(1.0, 3)
    d = m.^(-1/4)
    # d = fill(0.0, 3)
    N⁰ = fill(0.0, 3)
    K = size(σ,1)
    Λ = inv(σ)
    Q = (Λ + Λ')/2
    c = -Λ * d
    return EnergyConstraintProblem(σ,Λ,Q,c,m,d,N⁰,K,S,type)
end