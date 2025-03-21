""" Mainly for 2D visualization, help to understand how the geometry works and validate if the backend did the right thing and the interface correctly translate the problems"""
# """ Under construction """

function show_chevball(p::EnergyConstrProb, go_back::Bool=true)
    # !(p.K==2) && throw(ErrorException("can only show 2D"))
    itp = make_isotropic(p)
    if p.type == :indiv
        domain = InterPolySpheres(itp.A, itp.b, [], :indiv)
    elseif p.type == :total
        sp0 = Sphere(itp.yc, sqrt(itp.t))
        domain = InterPolySpheres(itp.A, itp.b, [sp0], :total)
    end
    sp = chevball(domain)
    
    color_sp = "orange"
    θ = 0:0.001:2π
    x = sp.c[1] .+ sp.r * cos.(θ)
    y = sp.c[2] .+ sp.r * sin.(θ)
    
    if go_back 
        color_sp = "cyan"
        x = itp.invL[1,1] * x .+ itp.invL[1,2] * y
        y = itp.invL[2,1] * x .+ itp.invL[2,2] * y
    end

    plot!(x, y, color=color_sp, linewidth=3)
end


function show_linear(
    p::Union{EnergyConstrProb,Nothing}=nothing,
    samples::Union{Vector{Vector{Float64}}, Nothing}=nothing,
)
    go_back = true # temporarily always plot in s-space
    go_back ? color_pts = "blue" : color_pts = "orange"

    # !isnothing(p) && plot!(create_poly(p), alpha=0.25, color="red")
    # !isnothing(p) && plot_sphere()

    !isnothing(p) && plot!(project(create_poly(p), [1,2]), alpha=0.35, color="purple", linewidth=3)
    !isnothing(samples) && scatter!(
        [s[1] for s in samples], [s[2] for s in samples], 
        color=color_pts, markersize=1.5, markerstrokewidth=0, ratio=1
    )
end

function show_quadratic(
    p::Union{EnergyConstrProb,Nothing}=nothing,
    samples::Union{Vector{Vector{Float64}}, Nothing}=nothing,
)
    go_back = true # temporarily always plot in s-space
    go_back ? color_pts = "blue" : color_pts = "orange"

    !isnothing(samples) && scatter!(
        [s[1] for s in samples], [s[2] for s in samples], 
        color=color_pts, markersize=1.5, markerstrokewidth=0, ratio=1
    )
end