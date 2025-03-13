"""
This script provides charaterization and visualization of energetic feasibility domains.
Currently it contains size of domain (Ω) and 2D visualization.
"""
#TODO: trophic profile (t), makie for fancy plots

using JuMP
using Polyhedra
using QHull
using Plots
using RecursiveArrayTools
include("_problem.jl")

"""Get the polyhedron directly from constraints (half-space inequalities)"""
function domain_by_constr(p::EnergyConstraintProblem)
    m, _, flag = create_model(p)
    if !flag
        throw(ArgumentError("direct constraint approach only applies to linear problem"))
    end
    return polyhedron(m)
end

"""Compute volume directly from constraints, only for linear problem."""
function vol_by_constr(p::EnergyConstraintProblem)
    poly = domain_by_constr(p)
    return volume(poly)
end


"""Compute volume from the convexhull of the sampling points."""
function vol_by_samp(samples::VectorOfArray)
    domain = vrep(Matrix(samples'))
    poly = polyhedron(domain, QHull.Library())
    removevredundancy!(poly)
    volume(poly)
end


"""Visualize the projection of EFD onto the i&j-th dimension"""
function visualize_2d(samples::VectorOfArray, p::EnergyConstraintProblem, dims::Vector{Int}=[1,2])
    if length(dims)!=2 || dims[1]>=dims[2]
        throw(ArgumentError("Incorrect visualization dimensions"))
    end
    
    # theoretical constraints
    m, _, flag = create_model(p)
    if(flag)
        poly_theory = polyhedron(m)
        proj_theory = project(poly_theory, dims)
        plot(proj_theory, alpha=0.5, color=:blue, label="constraints")
    end
    
    # convexhull
    domain = vrep(Matrix(samples'))
    poly_sampled = polyhedron(domain, QHull.Library())
    removevredundancy!(poly_sampled)
    proj_sampled = project(poly_sampled, dims)
    plot!(proj_sampled, alpha=0.6, color=:seagreen, label="convexhull")

    # sampling points
    x = Matrix(samples')[:, dims[1]]
    y = Matrix(samples')[:, dims[2]]
    scatter!(x, y, markersize=2, markerstrokewidth=0, alpha=0.5, color=:darkorange, label="samples")
    
    title!("Projection onto dimensions $dims")
    plot!(legend=:bottomright)
end