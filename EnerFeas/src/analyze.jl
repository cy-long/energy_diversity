# --- functional data analysis (after) --- 
function smooth(x::Vector{Float64}, y::Vector{Float64}, frac::Float64=0.02)
    if length(x) != length(y)
        throw(ArgumentError("x and y must have the same length"))
    end
    y[isinf.(y)] .= 0
    model = loess(x, y; span=frac)
    return model
end

function smooth_curve(x::Vector{Float64}, y::Vector{Float64}, frac::Float64=0.02; len::Int=1000)
    model = smooth(x, y, frac)
    xs = range(minimum(x), maximum(x), length=len)
    ys = predict(model, xs)
    return Vector(xs), ys
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
    y[isinf.(y)] .= 0.0
    x, y = smooth_curve(x, y, 0.02, len=length(x))

    initial = Int(floor(length(x) / 2));
    ind1 = (diff(y).>=0)
    ind2 = vcat(0,(diff(diff(y)) .<=0))
    ind2[1:initial] .= 1 # omit the first half

    score_0 = dot(abs.(diff(y)), ind1.*ind2) / sum(abs.(diff(y)))
    return score_0
end

function score_unimodal(x::Vector{Float64}, y::Vector{Float64})
    y[isinf.(y)] .= 0;

    x, y = smooth_curve(x, y, 0.02, len=length(x))

    xi, _ = extract_peak(x, y)
    dy = diff(y); xx = x[2:end]
    ind = ifelse.(dy .>= 0, xx .< xi, xx .> xi)
    score_1 = dot(abs.(dy), ind)/sum(abs.(dy))
    return score_1
end

# --- estimations of typical Q (before) ---

#> carefully test!!
function optimal_supply(p::EcosysParams)
    tp = translate_EFD(p, :matr) # include quadratic constraint
    A = p.N⁰' * (tp.cholP.U \ tp.yc);
    B = sqrt(p.N⁰' * (tp.P \ p.N⁰));
    C = 0.25 * tp.c' * (tp.P \ tp.c);
    B2 = B^2
    return A + B2/2 + sqrt((2*A + B2)^2 - 4*(A^2-B2*C))/2
end

function baseline_supply(p::EcosysParams)
    return dot(p.d + p.σ * (p.ϵ .* p.N⁰), p.N⁰)
end

function select_range(p::EcosysParams)
    Q0 = baseline_supply(p)
    Q1 = optimal_supply(p)
    Q_range = vcat(
        exp.(range(log(1e-4), log(Q0); length=25)),            # before critical
        exp.(range(log(Q0), log(Q1); length=351)[2:end]),     # around optimum
        exp.(range(log(Q1), log(max(10Q1, 100.0)); length=126)[2:end])  # tail
    )
    return Q_range
end