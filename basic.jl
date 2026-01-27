""" Basic tests for EnerFeas.jl """

using Pkg
Pkg.status()

using EnerFeas
using Random, Distributions, LinearAlgebra, Plots

#1 sampling initialization domain (a polyhedron)
p = generate_model_system(2, 24, 1.0, 1.25, 1.5; d_sd=0.1, σ0 = 0.125);
p.Q = 10.0;
check_feasible_EFD(p, :init, Bool.([true, true]))
tp = translate_EFD(p, :init);
samps = sample_EFD(tp, n_sample=10^4);
samps1 = sample_EFD(p, :init; n_sample=10^4);

plt_1 = plot(size=(300,300));
scatter!(plt_1, [s[1] for s in samps], [s[2] for s in samps], label="TransParams", markerstrokecolor=:auto);
scatter!(plt_1, [s[1] for s in samps1], [s[2] for s in samps1], label="Direct", alpha=0.5, color=:orange, markerstrokecolor=:auto);
display(plt_1)

#2 sampling quadratic constraint domain
check_feasible_EFD(p, :matr, Bool.([true, true]))
tp = translate_EFD(p, :matr);
samps = sample_EFD(tp, n_sample=10^4);

plt_2 = plot(size=(300,300));
scatter!(plt_2, [s[1] for s in samps], [s[2] for s in samps], label="TransParams", markerstrokecolor=:auto);
display(plt_2)


#3 volume estimations
volume_EFD(p, :init; exact=true)
volume_EFD(p, :matr)
Q_range = vcat(1.0:0.25:20.0);
vols = volume_range_EFD(p, :matr, Q_range)