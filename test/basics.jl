""" Testing basic functionality """

include("../src/EnerFeas.jl")
using .EnerFeas
using Random, Distributions, LinearAlgebra, Plots


#1 sampling linear constraint domain
ec1 = ecosys_config(K=3, S_type=:indiv, seed=24);
σ1 = generate_sigma_arrays(ec1, 1) 
p1 = generate_problem(ec1, σ1); p1.S = 4.5;

samples1 = sample_EFD(p1, 15000);

plot(ratio=1)
show_linear(p1, samples1)
show_chevball(p1)

#2 sampling quadratic constraint domain
ec2 = ecosys_config(K=2, S_type=:total, seed=12);
σ2 = generate_sigma_arrays(ec2, 1);
p2 = generate_problem(ec2, σ2); p2.S = 6.0;

samples2 = sample_EFD(p2, 15000);

plot(ratio=1)
show_quadratic(p2, samples2)

#3 linear volume estimation
volume_EFD(p1, true, 10, 1)
volume_EFD(p1, false, 10, 1)

#4 quadratic volume estimation
volume_EFD(p2, true, 10, 1) # error
volume_EFD(p2, false, 10, 1)