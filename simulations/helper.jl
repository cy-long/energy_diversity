""" Shared utilities for simulation scripts: CLI parsing, IO, math tools, and curve aggregation. """

using Statistics

# ── CLI & IO ──────────────────────────────────────────────────────────────────

function parse_arg(args::Vector{String}, idx::Int, default, ::Type{T}) where {T}
    return length(args) >= idx ? parse(T, args[idx]) : default
end

function unique_outfile(outdir::AbstractString, prefix::AbstractString, seed::Int)
    base = joinpath(outdir, "$(prefix)_seed$(seed)")
    outfile = base * ".jld2"
    k = 1
    while isfile(outfile)
        outfile = base * "_$(k).jld2"
        k += 1
    end
    return outfile
end

# ── Math utilities ────────────────────────────────────────────────────────────

logit01(x) = log(x / (1 - x))
expit(x) = inv(1 + exp(-x))

function logspace10(xmin::Real, xmax::Real, n::Int)
    return collect(10.0 .^ range(log10(xmin), log10(xmax), length=n))
end

function interp_linear_const(x::AbstractVector, y::AbstractVector, xs::AbstractVector)
    ys = similar(xs, Float64)
    j = 1
    n = length(x)

    for (i, xq) in pairs(xs)
        if xq <= x[1]
            ys[i] = y[1]
            continue
        end
        if xq >= x[end]
            ys[i] = y[end]
            continue
        end

        while j < n - 1 && x[j + 1] < xq
            j += 1
        end

        x0, x1 = x[j], x[j + 1]
        y0, y1 = y[j], y[j + 1]
        t = (xq - x0) / (x1 - x0)
        ys[i] = y0 + t * (y1 - y0)
    end

    return ys
end

# ── Curve aggregation ─────────────────────────────────────────────────────────

"""Aggregate probability-valued curves (e.g. Pₘ) via logit-space mean ± std."""
function estimate_functional_mean_std(
    rows;
    x_range=nothing,
    n_points::Int=300,
    epsilon::Float64=1e-3,
)
    isempty(rows) && return nothing

    if isnothing(x_range)
        xmin = minimum(first(row.Qs) for row in rows)
        xmax = maximum(last(row.Qs) for row in rows)
        x_range = (xmin, xmax)
    end

    xs = logspace10(x_range[1], x_range[2], n_points)
    y_interp = Vector{Vector{Float64}}()

    for row in rows
        length(row.Qs) < 2 && continue
        push!(y_interp, interp_linear_const(row.Qs, row.y, xs))
    end

    isempty(y_interp) && return nothing

    y_stack = reduce(vcat, permutedims.(y_interp))
    y_clamped = clamp.(y_stack, epsilon, 1.0 - epsilon)
    logit_y = logit01.(y_clamped)

    μ = mean(logit_y; dims=1) |> vec
    σ = std(logit_y; dims=1) |> vec

    return (
        xs=xs,
        mean=expit.(μ),
        lower=expit.(μ .- σ),
        upper=expit.(μ .+ σ),
        n_curves=size(y_stack, 1),
        x_range=x_range,
    )
end

"""Aggregate general curves (e.g. E[K|Q]) via ordinary mean ± std."""
function estimate_curve_mean_std(
    rows;
    x_range=nothing,
    n_points::Int=300,
)
    isempty(rows) && return nothing

    if isnothing(x_range)
        xmin = minimum(first(row.Qs) for row in rows)
        xmax = maximum(last(row.Qs) for row in rows)
        x_range = (xmin, xmax)
    end

    xs = logspace10(x_range[1], x_range[2], n_points)
    y_interp = Vector{Vector{Float64}}()

    for row in rows
        length(row.Qs) < 2 && continue
        push!(y_interp, interp_linear_const(row.Qs, row.y, xs))
    end

    isempty(y_interp) && return nothing

    y_stack = reduce(vcat, permutedims.(y_interp))
    μ = mean(y_stack; dims=1) |> vec
    σ = std(y_stack; dims=1) |> vec

    return (
        xs=xs,
        mean=μ,
        lower=max.(μ .- σ, 0.0),
        upper=μ .+ σ,
        n_curves=size(y_stack, 1),
        x_range=x_range,
    )
end
