using EnerFeas
using Plots
using LinearAlgebra

# include("../simulations/helper.jl")
include("../simulations/partition.jl")

function total_capture(s::AbstractVector{<:Real}, p)
    N = p.σ \ (Float64.(s) .- p.d)
    if N[1] > 1e-8 && N[2] > 1e-8
        return dot(Float64.(s), N)
    elseif N[1] > 1e-8
        return s[1] * (s[1] - p.d[1]) / p.σ[1, 1]
    elseif N[2] > 1e-8
        return s[2] * (s[2] - p.d[2]) / p.σ[2, 2]
    else
        return 0.0
    end
end

function boundary_rays(p)
    vs = Vector{Vector{Float64}}()
    for i in 1:2
        a, b = p.σ[i, i], p.d[i]
        λ = (sqrt(b^2 + 4a * p.Q) - b) / (2a)
        push!(vs, p.d + λ * p.σ[:, i])
    end
    return vs
end

p = generate_model_system(2, 100, 1.0, 1.2, 1.0)
p.Q = 20.0

plt = plot(
    ratio=1,
    xlim=(0, 8),
    ylim=(0, 8),
    xlabel="s₁",
    ylabel="s₂",
    title="maturation partitions (S=2)",
    legend=:topright,
    grid=false,
    colorbar=false,
)

comms, names = enumerate_comm(p.S; include_empty=false)
part_colors = Dict("12" => "#2e5f8a", "1" => "#e8a838", "2" => "#c25a4a")
alphas = create_alpha(names)

for (comm, name) in zip(comms, names)
    check_feasible_EFD(p, :matr, comm) || continue
    ss = sample_EFD(p, :matr, comm; n_sample=1500)
    scatter!(plt,
        [s[1] for s in ss],
        [s[2] for s in ss];
        color=get(part_colors, name, :gray90),
        alpha=alphas[name],
        label=name,
        markersize=3.5,
        markerstrokecolor=:auto,
        legend = :right
    )
end

# total capture = Q contour
x_grid = range(0, 10.0; length=200)
y_grid = range(0, 10.0; length=200)
z_grid = [total_capture([x, y], p) for y in y_grid, x in x_grid]
contour!(plt, x_grid, y_grid, z_grid;
    levels=[p.Q], color=:dodgerblue, linewidth=4, label=false, colorbar=false)

# dummy legend for energy boundary
plot!(plt, [], []; color=:dodgerblue, linewidth=4, label="Q = $(p.Q)")

# boundary rays and axes
for v in boundary_rays(p)
    plot!(plt, [p.d[1], v[1]], [p.d[2], v[2]];
        color=:gray50, linewidth=3, label=false)
end
plot!(plt, [0, p.d[1]], [p.d[2], p.d[2]]; color=:gray50, linewidth=3, label=false)
plot!(plt, [p.d[1], p.d[1]], [0, p.d[2]]; color=:gray50, linewidth=3, label=false)

display(plt)
savefig(plt, "figures/viz_partitions_2d.pdf")