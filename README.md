energy_diversity
====
This repository contains the Julia code for the paper: *Energetic Constraints determine the Diversity of Feasible Ecological Systems*

## Instructions
### 1. Install Julia
This project requires Julia ≥ 1.10. Download and install Julia from the [here](https://julialang.org/install/).

### 2. Install EnerFeas package
The scripts depend on the EnerFeas package, hosted in a separate [EnerFeas](https://github.com/cy-long/EnerFeas.jl) repository.
Activate the Julia environment at the root directory, then add the EnerFeas package:
```julia
using Pkg
Pkg.add(url="https://github.com/cy-long/EnerFeas.jl", rev="v0.3.0")
```
<!-- For ssh, use: Pkg.add(url="git@github.com:cy-long/EnerFeas.jl.git") -->

### 3. Set up other dependencies
In the Julia environment at the root directory, run the following command:
```julia
using Pkg
Pkg.activate(".")
Pkg.instantiate()
```
This will install all other necessary dependencies for the local environment.

### 4. Test basic functionality
The [basic.jl](https://github.com/cy-long/energy_diversity/blob/main/basic.jl) script provides a quick overview for define and compute feasibility domain with energetic constraint.

### 5. Reproduce figures and data
The dataset for theoretical analysis is stored in [data/output](https://github.com/cy-long/energy_diversity/tree/main/data/output), which can also be computed using [main/main.jl](https://github.com/cy-long/energy_diversity/blob/main/main/main.jl) and then exported using [main/theory.jl](https://github.com/cy-long/energy_diversity/blob/main/main/theory.jl).
Figure 2 is prepared under [main/fig2.py](https://github.com/cy-long/energy_diversity/blob/main/main/fig2.py);
Figure 3 is prepared under [main/fig3.py](https://github.com/cy-long/energy_diversity/blob/main/main/fig3.py);
All sensitivity analysis in the supplementary figures are prepared under scripts in [sensitivity](https://github.com/cy-long/energy_diversity/tree/main/sensitivity).
