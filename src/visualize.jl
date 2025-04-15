""" Mainly for 2D visualization, help to understand how the geometry works and validate if the backend did the right thing and the interface correctly translate the problems"""
# """ Under construction """

function show_chevball(p::EnergyConstrProb, go_back::Bool=true)
    # !(p.K==2) && throw(ErrorException("can only show 2D"))
    itp = translate_EFD(p)
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
    
    plt = plot()
    !isnothing(samples) && scatter!(
        plt,
        [s[1] for s in samples], [s[2] for s in samples], 
        color=color_pts, markersize=1.5, markerstrokewidth=0, ratio=1
    )
    return plt
end


function grow_quadratic(
    p0::Union{EnergyConstrProb,Nothing}=nothing,
    St::Float64=p0.S,
    plt=nothing
)
    if plt === nothing
        plt = plot()  # Create new plot if none provided
    end

    go_back = false
    color_line = go_back ? "blue" : "orange"

    p0.S = St

    L = cholesky(p0.Q).U; invL = inv(L)
    A = -vcat(I, p0.Λ) * invL
    b = -vcat(zeros(p0.K), p0.N⁰ - p0.c)
    yc = -0.5 * invL' * p0.c
    t = St + 0.25 * p0.c' * inv(p0.Q) * p0.c

    θ = 0:0.001:2π
    x_sp = yc[1] .+ sqrt(t) * cos.(θ)
    y_sp = yc[2] .+ sqrt(t) * sin.(θ)
    
    if go_back 
        x_sp_tmp = invL[1,1] * x_sp .+ invL[1,2] * y_sp
        y_sp = invL[2,1] * x_sp .+ invL[2,2] * y_sp
        x_sp = x_sp_tmp
    end

    plot!(plt, x_sp, y_sp, color=color_line, linewidth=3)

    samples = sample_EFD(p0, 1000, go_back)
    scatter!(
        plt,
        [s[1] for s in samples], [s[2] for s in samples], 
        color=color_line, markersize=1.5, markerstrokewidth=0, ratio=1
    )

    return plt
end