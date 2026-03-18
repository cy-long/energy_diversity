energy_diversity
====
This repository contains the Julia code for the paper: *Energetic Constraints determine the Diversity of Feasible Ecological Systems*

## Instructions
### 1. Install Julia
This project requires Julia ≥ 1.10. Download and install Julia from the [here](https://julialang.org/install/).

### 2. Set up dependencies
The scripts depend on the [EnerFeas](https://github.com/cy-long/EnerFeas.jl) package, which is not registered in Julia's default registry. On Julia 1.10, add `EnerFeas` from GitHub first, then instantiate the rest of the environment.

In the Julia environment at the root directory, run the following command:
```julia
using Pkg
Pkg.activate(".")
Pkg.add(url="https://github.com/cy-long/EnerFeas.jl", rev="v0.3.0")
Pkg.instantiate()
```
This will install the correct `EnerFeas` version together with the other dependencies for the local environment.

### 3. Test basic functionality
The [basic.jl](https://github.com/cy-long/energy_diversity/blob/main/basic.jl) script provides a quick overview for define and compute feasibility domain with energetic constraint.

### 4. Reproduce figures and data
The dataset for theoretical analysis is stored in [data/output](https://github.com/cy-long/energy_diversity/tree/main/data/output), which can also be computed using [main/main.jl](https://github.com/cy-long/energy_diversity/blob/main/main/main.jl) and then exported using [main/theory.jl](https://github.com/cy-long/energy_diversity/blob/main/main/theory.jl).
Figure 2 is prepared under [main/fig2.py](https://github.com/cy-long/energy_diversity/blob/main/main/fig2.py);
Figure 3 is prepared under [main/fig3.py](https://github.com/cy-long/energy_diversity/blob/main/main/fig3.py);
All sensitivity analysis in the supplementary figures are prepared under scripts in [sensitivity](https://github.com/cy-long/energy_diversity/tree/main/sensitivity).
