"""
Hello-world demo for EnerFeas.jl

For one ecosystem `p`, we vary total supply `Q` and 
compute the probability of maturation under chaging `Q`
"""

using EnerFeas
using Random
using Plots

# 1. Build a small reproducible ecosystem (seed=42).
p = generate_model_system(8, 42, 2.0, 1.0, 1.0)
comm = community(8, [1,3,4,6]) # {1,3,4,6}

# 2. Pick one reference Q and confirm the maturation domain exists.
reference_Q = 10.0
p.Q = reference_Q
@assert check_feasible_EFD(p, :matr)

# 3. Compute the maturation-domain volume at this reference Q.
volume_at_reference_Q = volume_EFD(p, :matr, comm)

# 4. Sweep over Q and compute the probability of maturation:
Q_range = select_range(p, comm);
matr_volumes = volume_range_EFD(p, :matr, Q_range, comm; n_sample=5*10^4, show_prog=true);
cons_volumes = volume_range_C(p, Q_range);
pm_curve = [c > 0 ? v / c : NaN for (v, c) in zip(matr_volumes, cons_volumes)];

ref_idx = findmin(abs.(Q_range .- reference_Q))[2];
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
);
scatter!(plt, [Q_range[ref_idx]], [pm_curve[ref_idx]]; ms=7, color=:red, label="ref. Q");
savefig(plt, "figures/basic_test_1.pdf")