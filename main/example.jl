""" A 2-sp example to illustrate the geometry of feasibility domains under supply gradient"""

using Pkg
Pkg.activate(".")

include("../src/EnerFeas.jl");
using .EnerFeas
using LinearAlgebra, LazySets, Optim
using CairoMakie, Plots, Colors
using JLD2

function customize_problem(ec::EcosysConfig, σ::Matrix{Float64})
    ϵ = [0.0,0.0];
    N⁰ = [1.0, 1.2];
    d = [1.2, 1.8];
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d;
    S = 0.0
    return EnergyConstrProb(σ,Λ,Q,c,ϵ,d,N⁰,ec.K,S,ec.S_type)
end;

function fill_region_with_dots!(region::HPolyhedron, ax::Axis, a::Float64)
    dx, dy = 0.27, 0.225
    X = 0.0:dx:1.2a
    Y = 0.0:dy:a
    grid_points = [(x,y) for x in X, y in Y]

    points = vec(grid_points)
    inside = [p for p in points if collect(p) ∈ region]

    Xs = [p[1] for p in inside]
    Ys = [p[2] for p in inside]

    CairoMakie.scatter!(ax, Xs, Ys; color =(:darkorange,1.0), markersize = 1.55)
end;

function visualize_ellipsoid!(p::EnergyConstrProb, ax::Axis)
    ax.backgroundcolor = :transparent
    itp = translate_EFD(p)
    QQ = (itp.t * inv(p.Q) + (itp.t * inv(p.Q))')/2
    E = Ellipsoid(itp.invL * itp.yc, QQ)
    P = overapproximate(E, PolarDirections(200));
    verts = vertices_list(P);
    X = [v[1] for v in verts];
    Y = [v[2] for v in verts];
    push!(X, X[1]); push!(Y, Y[1]);
    poly!(
        ax,
        X, Y;
        color = (:white, 0),       # fully transparent fill
        strokecolor = :dodgerblue, # boundary color
        strokewidth = 1            # boundary thickness
)
end;

function visualize_cone!(p::EnergyConstrProb, a::Float64, ax::Axis)
    λ1 = minimum(filter(>(0), vcat(([1.2a,a].-p.d)./p.σ[:,1], (-p.d)./p.σ[:,1])))
    λ2 = minimum(filter(>(0), vcat(([1.2a,a].-p.d)./p.σ[:,2], (-p.d)./p.σ[:,2])))
    lines!(ax, [p.d[1], p.d[1]+λ1*p.σ[1,1]], [p.d[2], p.d[2]+λ1*p.σ[2,1]], color = :gray50, linewidth = 1)
    lines!(ax, [p.d[1], p.d[1]+λ2*p.σ[1,2]], [p.d[2], p.d[2]+λ2*p.σ[2,2]], color = :gray50, linewidth = 1)
    CairoMakie.scatter!(ax, p.d..., color=:gray50, markersize=5);
end;

function visualize_linear!(p::EnergyConstrProb, ax::Axis)
    p1 = [p.S / p.N⁰[1], 0.0]
    p2 = [0.0, p.S / p.N⁰[2]]
    lines!(ax, [p1[1], p2[1]], [p1[2], p2[2]], color = :darkorange, linewidth = 1)
end;

function visualize_feasdom!(p::EnergyConstrProb, ax::Axis, a::Float64, cl=:gray70)
    itp = translate_EFD(p)
    QQ = (itp.t * inv(p.Q) + (itp.t * inv(p.Q))')/2
    E = Ellipsoid(itp.invL * itp.yc, QQ);
    A = -vcat(I, p.Λ, -p.N⁰')
    b = -vcat(zeros(p.K), p.ϵ .* p.N⁰ - p.c, -p.S)
    PD = HalfSpace([1.0, 0.0], 1.2a) ∩
         HalfSpace([-1.0, 0.0], 0.0) ∩
         HalfSpace([0.0, 1.0], a) ∩
         HalfSpace([0.0, -1.0], 0.0) ∩
         HalfSpace(p.N⁰, p.S)
    FI = HPolyhedron(A, b);
    FR = E ∩ FI;

    P = overapproximate(FR, PolarDirections(1000))
    verts = vertices_list(P)
    X = [v[1] for v in verts]
    Y = [v[2] for v in verts]
    push!(X, X[1]); push!(Y, Y[1])
    # visualize PD here with gray90
    P_PD = overapproximate(PD, PolarDirections(500))
    verts_PD = vertices_list(P_PD)
    X_PD = [v[1] for v in verts_PD]
    Y_PD = [v[2] for v in verts_PD]
    push!(X_PD, X_PD[1]); push!(Y_PD, Y_PD[1])
    poly!(ax, X_PD, Y_PD; color=cl, strokecolor=(:white, 0), strokewidth=1)
    poly!(ax, X, Y; color=(:dodgerblue, 0.5), strokecolor=(:white, 0), strokewidth=1)
    fill_region_with_dots!(FI, ax, a)
end;

function visualize_EFD(p::EnergyConstrProb, a::Float64, cl=:gray70)
    fig = Figure(size=(150,150), backgroundcolor=:transparent)
    ax = Axis(fig[1,1];
        aspect = DataAspect(),
        rightspinevisible = false,
        topspinevisible = false,
        xgridvisible = false,
        ygridvisible = false,
        xticksvisible = false,
        yticksvisible = false,
        xticklabelsvisible = false,
        yticklabelsvisible = false,
        xlabel = "",
        ylabel = "",
        xlabelsize = 8,
        ylabelsize = 8,
        xlabelpadding = 5,
        ylabelpadding = 5,
        backgroundcolor = :transparent,
    )
    CairoMakie.xlims!(ax, 0, 1.2a);
    CairoMakie.ylims!(ax, 0, a);
    visualize_ellipsoid!(p, ax);
    visualize_cone!(p, a, ax);
    visualize_linear!(p, ax);
    visualize_feasdom!(p, ax, a, cl);

    return fig
end;

# --- create a 2-sp example with seed=35 ---
ec = ecosys_config(K=2, S_type=:total, seed=35, n_scale=2.0, positive=true);
σ = generate_sigma_arrays(ec, 1);
p = customize_problem(ec, σ);
a = 10.0; # canvas size

p.S = 15.0
pl = visualize_EFD(p, a, ());
CairoMakie.save("figures/geom_Q3_1.svg",pl)


# --- calculate probabilities ---
p = customize_problem(ec, σ);
p1 = customize_problem(ec, σ); p1.type=:indiv;
Q_range = vcat(1e-5:0.1:20.0,20:0.5:100.0,100.0:5.0:1000.0);
vols = volume_range_EFD(p, Q_range, n_sample=3*10^4); # @load "data/example.jld2"
devols = [Q^2/(2*prod(p.N⁰)) for Q in Q_range];

vols1 = volume_range_EFD(p1, Q_range); 
Qc, _ = critical_energy(p);
Qopt,_ = extract_peak(Q_range, vols./devols, 0.02);

# @save "data/example.jld2" Q_range vols vols1 devols p p1 Qc Qopt;

# --- plot the results ---
plt1 = Plots.plot(framestyle=:box, grid=false,
    xlabel="",
    ylabel="",
    xaxis = :log10,
    xlim=(1,250),
    ylim = (0.0, 0.32),
    yticks = [0.1, 0.2],
    guidefont=font("Helvetica", 6),
    tickfont=font("Helvetica", 5),
    legendfont=font("Helvetica", 5),
    legend=:right,
    size=(240, 90),
    lw=1.0,
    foreground_color_legend = nothing,
    legend_border = false,
    background_color = :transparent,
);

Plots.plot!(plt1, Q_range, vols1./devols, color=:darkorange, label="", linewidth=1.5, linestyle=:dot);
Plots.plot!(plt1, Q_range, vols./devols, color=:dodgerblue, label="",linewidth=1.5,linestyle=:solid);

for (Qd,cl) in zip([7,15], [:gray90, :gray20])
    Plots.vline!(plt1, [Qd], color=cl, linewidth=1.0, label=false, alpha=0.2);
end

display(plt1)
save("figures/single.svg", plt1);