""" (Draft) A 2D example of increasing Q, visualize the probability functions and the geometry of domains """

using Revise
includet("../src/EnerFeas.jl");
using .EnerFeas
using LazySets, CairoMakie, LinearAlgebra
using Plots
gr();  # ensure GR backend
default(fontfamily = "Helvetica")
using Optim
using JLD2

function customize_problem(ec::EcosysConfig, σ::Matrix{Float64})
    ϵ = [0.0,0.0];
    N⁰ = [1.0, 1.35];
    d = [1.2, 1.8];
    Λ = inv(σ); Q = (Λ + Λ')/2; c = -Λ * d;
    S = 0.0
    return EnergyConstrProb(σ,Λ,Q,c,ϵ,d,N⁰,ec.K,S,ec.S_type)
end;

# create a real problem
ec = ecosys_config(K=2, S_type=:total, seed=35, n_scale=2.0, positive=true);
σ = generate_sigma_arrays(ec, 1);
p = customize_problem(ec, σ);
a = 11.0; # canvas size

function fill_region_with_dots!(region::HPolyhedron, ax::Axis, a::Float64)
    dx, dy = 0.27, 0.225
    X = 0.0:dx:(1.2*a)
    Y = 0.0:dy:a
    grid_points = [(x,y) for x in X, y in Y]

    points = vec(grid_points)
    inside = [p for p in points if collect(p) ∈ region]

    Xs = [p[1] for p in inside]
    Ys = [p[2] for p in inside]

    CairoMakie.scatter!(ax, Xs, Ys; color = (:darkorange, 0.9), markersize = 1.45)
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
    λ1 = minimum(filter(>(0), vcat(([1.2*a,a].-p.d)./p.σ[:,1], (-p.d)./p.σ[:,1])))
    λ2 = minimum(filter(>(0), vcat(([1.2*a,a].-p.d)./p.σ[:,2], (-p.d)./p.σ[:,2])))
    lines!(ax, [p.d[1], p.d[1]+λ1*p.σ[1,1]], [p.d[2], p.d[2]+λ1*p.σ[2,1]], color = :gray50, linewidth = 1)
    lines!(ax, [p.d[1], p.d[1]+λ2*p.σ[1,2]], [p.d[2], p.d[2]+λ2*p.σ[2,2]], color = :gray50, linewidth = 1)
    CairoMakie.scatter!(ax, p.d..., color=:gray50, markersize=5);
end;

function visualize_linear!(p::EnergyConstrProb, ax::Axis)
    p1 = [p.S / p.N⁰[1], 0.0]
    p2 = [0.0, p.S / p.N⁰[2]]
    lines!(ax, [p1[1], p2[1]], [p1[2], p2[2]], color = :darkorange, linewidth = 1)
end;

function visualize_feasdom!(p::EnergyConstrProb, ax::Axis, a::Float64)
    itp = translate_EFD(p)
    QQ = (itp.t * inv(p.Q) + (itp.t * inv(p.Q))')/2
    E = Ellipsoid(itp.invL * itp.yc, QQ);
    A = -vcat(I, p.Λ, -p.N⁰')
    b = -vcat(zeros(p.K), p.ϵ .* p.N⁰ - p.c, -p.S)
    FI = HPolyhedron(A, b);
    FR = E ∩ FI;

    P = overapproximate(FR, PolarDirections(1000))
    verts = vertices_list(P)
    X = [v[1] for v in verts]
    Y = [v[2] for v in verts]
    push!(X, X[1]); push!(Y, Y[1])
    poly!(ax, X, Y; color=(:dodgerblue, 0.8), strokecolor=(:white,0), strokewidth=1)
    fill_region_with_dots!(FI, ax, a)
end;

function visualize_EFD(p::EnergyConstrProb, a::Float64, i::Int)
    x_alphas = [0.0, 0.0, 0.0];
    y_alphas = [0.0, 0.0, 0.0];
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
        xlabel = "Supply 1",
        ylabel = "Supply 2",
        xlabelcolor=(:black, x_alphas[i]),
        ylabelcolor=(:black, y_alphas[i]),
        xlabelsize = 10,
        ylabelsize = 10,
        xlabelpadding = 5,
        ylabelpadding = 5,
        backgroundcolor = :transparent
    )
    CairoMakie.xlims!(ax, 0, 1.2*a);
    CairoMakie.ylims!(ax, 0, a);
    visualize_ellipsoid!(p, ax);
    visualize_cone!(p, a, ax);
    visualize_linear!(p, ax);
    visualize_feasdom!(p, ax, a);

    return fig
end;

Q_values = [5.0, 11.011021, 20.0];
for (i, Q) in enumerate(Q_values)
    p.S = Q;
    fig = visualize_EFD(p, a, i);
    display(fig);
    CairoMakie.save("figures/geometry_Q_$(i).pdf", fig);
end

# -------- Create a dummy legend --------
plt3 = Plots.plot(
    grid=false, framestyle=:none, legend=:bottomleft, legendfont=font("Helvetica", 4), background_color=:transparent, foreground_color_legend = nothing, legend_border = false, size=(150,50)
);

Plots.scatter!(plt3, [0.0], [0.0]; color=:white, alpha=0, markerstrokecolor=:auto, markersize=0.1, label="");
Plots.scatter!(plt3, [NaN], [NaN]; color=:darkorange, marker=:circle, markerstrokecolor=:auto, markersize=0.25, label="Initiation Domain");
Plots.scatter!(plt3, [NaN], [NaN]; color=:dodgerblue, marker=:rect, markerstrokecolor=:auto, markersize=1, strokealpha = 0,label="Realization Domain");

display(plt3)

savefig(plt3, "figures/dummy.pdf");

# -------- Generate the volume curves --------
Q_range = vcat(1e-5:0.1:20.0,20:0.5:100.0,100.0:5.0:1000.0);
# vols = volume_range_EFD(p, Q_range, n_sample=3*10^4);
devols = [Q^2/(2*prod(p.N⁰)) for Q in Q_range];

# p1 = customize_problem(ec, σ); p1.type=:indiv;
# vols1 = volume_range_EFD(p1, Q_range); 
@load "data/single.jld";
Qc, _ = critical_energy(p);
Qopt,_ = extract_peak(Q_range, vols./devols, 0.02);


plt1 = Plots.plot(framestyle=:box, grid=false,
    xlabel="Upper Energy Bound ",
    ylabel="Probability",
    xaxis = :log10,
    xlim=(1,250),
    ylim = (0.0, 0.32),
    yticks = [0.15, 0.3],
    guidefont=font("Helvetica", 6),
    tickfont=font("Helvetica", 5),
    legendfont=font("Helvetica", 5),
    legend=(0.80,0.40),
    size=(380, 120),
    lw=1.0,
    foreground_color_legend = nothing,
    legend_border = false,
    background_color = :transparent,
);

Plots.vspan!(plt1, [1.0, Qc], color=:white, alpha=0.1 ,label=false);
Plots.vspan!(plt1, [Qc, Qopt], color=:peachpuff, alpha=0.35, label=false);
Plots.vspan!(plt1, [Qopt, 250.0], color=:skyblue, alpha=0.25,label=false);

Plots.vline!(plt1, [Qc], color=:gray70, linewidth=0.8, linestyle = :solid, label=false);
Plots.vline!(plt1, [Qopt], color=:gray70, linewidth=0.8, linestyle = :solid, label=false);

Plots.plot!(plt1, Q_range, vols1./devols, color=:darkorange, label="Prob. Initiation", linewidth=1.5);
Plots.plot!(plt1, Q_range, vols./devols, color=:dodgerblue, label="Prob. Realization",linewidth=1.5);

for Qd in Q_values
    Plots.vline!(plt1, [Qd], color=:gray40, linewidth=0.8, linestyle = :dash, label=false);
end

display(plt1)

save("figures/single.pdf", plt1);


""" Some thoughts as we are getting clear about the picture """
# 1. The biomass accumulation process is a natural way to introduce all the conditions/constraints. N0, minimally viable biomass, is minimal in a sense that can "start the biomass accumulation" (instrinic growthrate) -> Methods section

# 2. Before jumping into calculating the probability, build some intuition and understanding. For instance, which of the constraint is effectively constraining the system (initially N0, then N*); initially the s∈D^S cannot generate a positive steady state; all that s inside the cone (N^*>0) takes much more energy to support a minimal system; as long as Q>0, one can arbitrarily create a feasible steady state, but the N^* should be arbitrarily small... These pictures and processes provide intuition in ecological terms, build intuition for readers. Probability is one formal way to aggregate these results -> Results section

## 2.1 The geometry is how we approach these pictures, therefore adding these thinkings into the main Figures (1 or 2). They are self-explaining, but we need to inteprete. 

# 3. Perhaps one related analysis: how does <N*> change when Q increases? How does that compare to N0? -> Consistency

# 4. About the critical Q thing: in our framework, N⁰s≤Q is what makes this nonzero critical Q; in fact ϵ=0 doesn't make a difference. However, if removing the condition of N⁰s≤Q in Feasibility domain, then any small Q lead to non-empty domains (we are quite clear about the math now: since d'Pd+c'd=0), adding ϵ>0 provides a nonzero basline to Qc, and perhaps certain σ could lower Qc to this baseline (solving Serguei's puzzle)  

# 5. Further idea: if 1-3 makes perfect since, how can we couple these energetic analysis and the biomass dynamics? How is this related to those theory papers that says available resources could structure the system? Considering they find phases only in terms of N*, can our framework provide more subtlty? But I don't know how to couple this changingQ-Ndynamics. Get really irrelavant.