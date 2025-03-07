include("src/lvmodel.jl")
include("src/sampler.jl")
include("src/domain.jl")
const GRB_ENV = Gurobi.Env(output_flag=0)

# Define the energetic problem
p = create_energy_problem(5, 8.0, :Individual, seed=20);
baseline_supply(p)

# Volume by triangularization
vol_by_constr(p)

# Volume by sampling and convexhull --> Currently not very efficient... A better way?
samp = create_sampler(p);
warmup!(samp)
samples = hr_sample!(samp, 20000);
vol_by_samp(samples)

# Visulization
visualize_2d(samples, p, [2,3])

# Energetic Trajectory of sampling
tot_supp = [total_supply(sample, p) for sample in samples]
ind_supp = [individual_supply(sample, p) for sample in samples]

plot(range(1, length=length(tot_supp)), tot_supp, xlabel="sampling run", ylabel="sampling supply", label="S")