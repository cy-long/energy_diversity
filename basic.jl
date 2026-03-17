"""
Hello-world demo for EnerFeas.jl

For one ecosystem `p`, we vary total supply `Q` and 
compute the probability of maturation under chaging `Q`
"""

using EnerFeas
using Random
using Plots

# 1. Build a small reproducible ecosystem (seed=42).
p = generate_model_system(3, 42, 2.0, 1.0, 1.0)

# 2. Pick one reference Q and confirm the maturation domain exists.
reference_Q = 10.0
p.Q = reference_Q
@assert check_feasible_EFD(p, :matr)

# 3. Compute the maturation-domain volume at this reference Q.
volume_at_reference_Q = volume_EFD(p, :matr)

# 4. Sweep over Q and compute the probability of maturation:
Q_range = collect(10 .^ range(-0.1, 2.0, length=150))
matr_volumes = volume_range_EFD(p, :matr, Q_range; n_sample=5*10^4, show_prog=false)
constraint_volumes = volume_range_C(p, Q_range)
pm_curve = [c > 0 ? v / c : NaN for (v, c) in zip(matr_volumes, constraint_volumes)]
reference_idx = findmin(abs.(Q_range .- reference_Q))[2]
p.Q = reference_Q

plt = plot(
    Q_range,
    pm_curve;
    xaxis=:log10,
    xlabel="Total supply Q",
    ylabel="Probability of maturation",
    title="PM-Q trajectory for a single ecosystem",
    lw=2.5,
    color=:blue,
    label="prob. maturation",
    size=(640, 400),
)
scatter!(plt, [Q_range[reference_idx]], [pm_curve[reference_idx]]; ms=7, color=:red, label="reference Q")
display(plt)