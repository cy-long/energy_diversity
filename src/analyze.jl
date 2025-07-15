using LinearAlgebra, Loess

# tackle the problem of inf in y

function smooth(x::Vector{Float64}, y::Vector{Float64}, frac::Float64=0.02)
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    y[isinf.(y)] .= 0
    model = loess(x, y; span=frac)
    return model
end

function smooth_curve(x::Vector{Float64}, y::Vector{Float64}, frac::Float64=0.02)
    model = smooth(x, y, frac)
    xs = range(minimum(x), maximum(x), length=1000)
    ys = predict(model, xs)
    return xs, ys
end

#> may also be extracted before the characterization: just an optimization problem
function extract_critical(x::Vector{Float64}, y::Vector{Float64})
    if y[end] <= 1e-7
        return NaN
    end
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    y[isinf.(y)] .= 0
    x_c = x[findall(y.==0)[end]]
    return x_c
end

function extract_peak(x::Vector{Float64}, y::Vector{Float64}, frac=0.02)
    if abs(maximum(y)) <= 1e-7
        return NaN
    end
    xs, ys = smooth_curve(x, y, frac)
    return xs[argmax(ys)], maximum(ys)
end

function score_saturation(x::Vector{Float64}, y::Vector{Float64})
    y[isinf.(y)] .= 0
    ind = (diff(y).>=0);
    score_0 = dot(abs.(diff(y)), ind)/sum(abs.(diff(y)))
    return score_0
end

function score_unimodal(x::Vector{Float64}, y::Vector{Float64})
    y[isinf.(y)] .= 0;
    xi, _ = extract_peak(x, y)
    dy = diff(y); xx = x[2:end]
    ind = ifelse.(dy .>= 0, xx .< xi, xx .> xi)
    score_1 = dot(abs.(dy), ind)/sum(abs.(dy))
    return score_1
end

function critical_energy(p::EnergyConstrProb)
    model = Model()
    @variable(model, s[i = 1:p.K]>=0)
    @variable(model, Qc >= 0)

    A = -vcat(I, p.Λ, -p.N⁰')
    b = -vcat(zeros(p.K), p.ϵ .* p.N⁰ - p.c, -Qc)
    @constraint(model, A*s .<= b)

    if p.type == :total
        @constraint(model, s'*p.Q*s + p.c'*s <= Qc)
    end

    @objective(model, Min, Qc)
    set_optimizer(model, SCS.Optimizer);
    set_optimizer_attribute(model, "verbose", 0)
    optimize!(model)
    if termination_status(model) != MOI.OPTIMAL
        return false
    else 
        return objective_value(model), value.(s)
    end
end

function optimal_supply(p::EnergyConstrProb)
    itp = translate_EFD(p)
    A = p.N⁰' * itp.invL * itp.yc;
    B = sqrt(p.N⁰' *inv(p.Q)* p.N⁰);
    C = 0.25 * p.c' * inv(p.Q) * p.c;
    B2 = B^2
    return A + B2/2 + sqrt((2*A + B2)^2 - 4*(A^2-B2*C))/2
end

# Calculate the baseline supply needed to sustain the ecosystem with least biomass
baseline_supply(p::EnergyConstrProb) = dot(p.d + p.σ * (p.ϵ .* p.N⁰), p.N⁰)

# Calculate the total supply at state s (weighed by N)
total_supply(s::Vector{Float64}, p::EnergyConstrProb) = transpose(s) * p.Λ * (s - p.d)