""" Profiling benchmark for EFD volume computation. """

using Revise
using Profile
using StatProfilerHTML
using EnerFeas

S = 8
seed = 42
p = generate_model_system(S, seed, 1.0, 1.0, 1.0);
Q_range = collect(10 .^ range(-0.1, 3.1, length=500));

vols_c = volume_range_C(p, Q_range);
vols_m = volume_range_EFD(p, :matr, Q_range);


vols_i1 = volume_range_EFD(p, :init, Q_range; exact=false);

plot(Q_range, vols_i./vols_c, label="init", xscale=:log10)
plot!(Q_range, vols_i1./vols_c, label="init1", xscale=:log10)

#= Profile.clear()
@profile begin
    res = build_partition_geometry_family(p; Q_range=Q_range,prog_dt=0.1);
end =#

statprofilehtml()