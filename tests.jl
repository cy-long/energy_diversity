""" testing """

include("src/EnerFeas.jl")
using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots


m0 = fill(1.0,3)

#1 sampling linear constraint domain
σ1 = generate_trophic("chain")
p1 = test_linear(σ1, m0, 3.5);
samples1 = sample_EFD(p1, 15000);

plot(ratio=1)
show_linear(p1, samples1)
show_chevball(p1)

#2 sampling quadratic constraint domain
σ2 = generate_trophic("exploitative")
p2 = test_quadratic(σ2, m0, 6.0, 0.75);
samples2 = sample_EFD(p2, 15000);

plot(ratio=1)
show_quadratic(p2, samples2)

#3 linear volume estimation
volume_EFD(p1, true, 10, 1)
volume_EFD(p1, false, 10, 1)

#4 quadratic volume estimation
volume_EFD(p2, true, 10, 1) # error
volume_EFD(p2, false, 10, 1)