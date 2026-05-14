[![DOI](https://zenodo.org/badge/940171603.svg)](https://doi.org/10.5281/zenodo.20185062)

energy_diversity
====
Code for the paper: *Energetic Constraints Shape the Diversity of Feasible Ecological Networks*

## Setup

**Julia ≥ 1.10** required. Install from [julialang.org](https://julialang.org/install/).

The core dependency [EnerFeas.jl](https://github.com/cy-long/EnerFeas.jl) is not in the default registry. Set up the environment with:
```julia
using Pkg
Pkg.activate(".")
Pkg.add(url="https://github.com/cy-long/EnerFeas.jl", rev="v0.3.1")
Pkg.instantiate()
```

Run `basic.jl` to verify the installation — it computes a single PM–Q curve and saves a test figure.

## Repository structure

```
simulations/       HPC numerical simulations (Julia)
  fig3/              run.jl / collect.jl — feasibility domains (Fig. 3)
  fig4/              run.jl / collect.jl — feasibility partitions (Fig. 4)
  helper.jl          helper functions
  partition.jl       partition-specific functions

plots/             Final figure generation (Python)
  fig2.py, fig3.py, fig4.py

sensitivity/       Supplementary sensitivity analyses (Julia, self-contained)

data/output/       Pre-computed datasets ready for plotting
```

## Running simulations

Each `simulations/fig*/` folder contains a `run.jl` (HPC computation) and a `collect.jl` (aggregate results → JSON for plotting).

**Fig. 3** — feasibility domain (full community) and parameter scalings:
```bash
julia --project=. simulations/fig3/run.jl <seed>
julia --project=. simulations/fig3/collect.jl
```

**Fig. 4** — feasibility partitions (all candidate communities):
```bash
julia --project=. simulations/fig4/run.jl <seed> [σsc] [d0] [N0] [n_chains] [n_sample] [n_layer] [show_prog] [outdir]
julia --project=. simulations/fig4/collect.jl
```

After `collect.jl` exports JSON to `data/output/`, use the corresponding Python scripts in `plots/` to generate final figures.
