""" Collect partition results across seeds and export data for Fig. 4 plotting. """

using Statistics
using DataFrames, Glob, JLD2
using JSON
using Plots

include(joinpath(@__DIR__, "..", "helper.jl"))

function collect_all_results(dir::String)
    files = sort(
        glob("partition_seed*.jld2", dir);
        by=file -> parse(Int, match(r"partition_seed(\d+)\.jld2$", basename(file)).captures[1]),
    )
    system_tables = DataFrame[]
    community_tables = DataFrame[]

    for file in files
        @info "Loading $file"
        push!(system_tables, load(file, "df_system"))
        push!(community_tables, load(file, "df_community"))
    end

    return (
        files=files,
        df_system=vcat(system_tables...),
        df_community=vcat(community_tables...),
    )
end

function community_colors(keys::Vector{String})
    palette = [
        "#6C8EBF", # muted blue
        "#C97B84", # dusty rose
        "#7BAE7F", # sage green
        "#B48EAD", # muted purple
        "#C2A35D", # soft ochre
        "#6FA3A6", # muted teal
        "#D08770", # muted salmon
        "#8F9779", # olive gray
        "#7C9FB0", # blue gray
        "#A89078", # taupe
        "#9C89B8", # lavender gray
        "#5E8C6A", # moss green
    ]
    colors = Dict(key => palette[mod1(i, length(palette))] for (i, key) in pairs(keys))
    if "∅" in keys
        colors["∅"] = "#D9D9D9"
    end
    return colors
end

function ordered_comm_keys(df_community::DataFrame; S::Union{Nothing, Int}=nothing)
    table = select(df_community, :S, :comm_key, :comm_size) |> unique
    if !isnothing(S)
        table = subset(table, :S => ByRow(isequal(S)))
    end
    sort!(table, [:comm_size, :comm_key])
    return collect(table.comm_key)
end

function build_pm_curve_table(df_system::DataFrame, df_community::DataFrame)
    curves = leftjoin(
        select(df_community, :seed, :S, :comm_key, :comm_size, :pm),
        select(df_system, :seed, :S, :Q_range),
        on=[:seed, :S],
    )
    rename!(curves, :Q_range => :Qs, :pm => :y)
    return curves
end

function plot_pm_functional_by_size(
    df_system::DataFrame,
    df_community::DataFrame;
    S::Int=4,
    sizes::Union{Nothing, Vector{Int}}=nothing,
    n_points::Int=300,
    epsilon::Float64=1e-3,
    xmin::Float64=1.0,
    xmax::Float64=150.0,
    ymax::Float64=0.35,
    band_alpha::Float64=0.18,
    linewidth::Float64=2.0,
    plot_size::Tuple{Int, Int}=(520, 360),
)
    curves = build_pm_curve_table(df_system, df_community)
    curves = subset(curves, :S => ByRow(isequal(S)))

    if isnothing(sizes)
        sizes = sort(unique(subset(df_community, :S => ByRow(isequal(S))).comm_size))
    end

    size_colors = community_colors(string.(sizes))
    plt = plot(
        size=plot_size,
        legend=:topright,
        grid=false,
        xaxis=:log10,
        xlims=(xmin, xmax),
        ylims=(0.0, ymax),
        xlabel="Q",
        ylabel="Pₘ",
        title="Average Pₘ by size; S=$(S)",
        framestyle=:box,
        legend_foreground_color=nothing,
    )

    for comm_size in sizes
        sub = subset(curves, :comm_size => ByRow(isequal(comm_size)))
        rows = collect(eachrow(sub))
        stats = estimate_functional_mean_std(
            rows;
            x_range=(xmin, xmax),
            n_points=n_points,
            epsilon=epsilon,
        )
        isnothing(stats) && continue

        color = size_colors[string(comm_size)]
        plot!(
            plt,
            stats.xs,
            stats.lower;
            fillrange=stats.upper,
            fillalpha=band_alpha,
            linewidth=0,
            color=color,
            label=nothing,
        )
        plot!(
            plt,
            stats.xs,
            stats.mean;
            color=color,
            linewidth=linewidth,
            label="size $(comm_size) (n=$(stats.n_curves))",
        )
    end

    return plt
end

function plot_pm_functional_by_size_weighted(
    df_system::DataFrame,
    df_community::DataFrame;
    S::Int=4,
    sizes::Union{Nothing, Vector{Int}}=nothing,
    n_points::Int=300,
    epsilon::Float64=1e-3,
    xmin::Float64=1.0,
    ymax=nothing,
    band_alpha::Float64=0.18,
    linewidth::Float64=2.0,
    plot_size::Tuple{Int, Int}=(520, 360),
)
    curves = build_pm_curve_table(df_system, df_community)
    curves = subset(curves, :S => ByRow(isequal(S)))

    if isnothing(sizes)
        sizes = sort(unique(subset(df_community, :S => ByRow(isequal(S))).comm_size))
    end

    xmax = maximum(last(q) for q in curves.Qs)
    size_colors = community_colors(string.(sizes))
    ymax_auto = 0.0
    weighted_stats = NamedTuple[]

    for comm_size in sizes
        sub = subset(curves, :comm_size => ByRow(isequal(comm_size)))
        rows = collect(eachrow(sub))
        n_types = length(unique(sub.comm_key))
        stats = estimate_functional_mean_std(
            rows;
            x_range=(xmin, xmax),
            n_points=n_points,
            epsilon=epsilon,
        )
        isnothing(stats) && continue

        weighted = (
            comm_size=comm_size,
            n_types=n_types,
            xs=stats.xs,
            mean=n_types .* stats.mean,
            lower=n_types .* stats.lower,
            upper=n_types .* stats.upper,
        )
        ymax_auto = max(ymax_auto, maximum(weighted.upper))
        push!(weighted_stats, weighted)
    end

    ymax_val = isnothing(ymax) ? ymax_auto : ymax
    plt = plot(
        size=plot_size,
        legend=:topright,
        grid=false,
        xaxis=:log10,
        xlims=(xmin, xmax),
        ylims=(0.0, ymax_val),
        xlabel="Q",
        ylabel="n × Pₘ",
        title="Weighted average Pₘ; S=$(S), agg. by size",
        framestyle=:box,
        legend_foreground_color=nothing,
    )

    for stats in weighted_stats
        color = size_colors[string(stats.comm_size)]
        plot!(
            plt,
            stats.xs,
            stats.lower;
            fillrange=stats.upper,
            fillalpha=band_alpha,
            linewidth=0,
            color=color,
            label=nothing,
        )
        plot!(
            plt,
            stats.xs,
            stats.mean;
            color=color,
            linewidth=linewidth,
            label="size $(stats.comm_size) (n=$(stats.n_types))",
        )
    end

    return plt
end

function compute_expected_richness_curve(
    df_system::DataFrame,
    df_community::DataFrame;
    S::Int,
    seed::Int,
)
    sys = subset(df_system, :S => ByRow(isequal(S)), :seed => ByRow(isequal(seed)))
    nrow(sys) == 1 || error("Expected exactly one system row for S=$(S), seed=$(seed), got $(nrow(sys))")

    comm = subset(df_community, :S => ByRow(isequal(S)), :seed => ByRow(isequal(seed)))
    Qs = sys.Q_range[1]
    EK = zeros(length(Qs))

    for row in eachrow(comm)
        EK .+= row.comm_size .* row.pm
    end

    return (
        seed=seed,
        S=S,
        Qs=Qs,
        y=EK,
        n_communities=nrow(comm),
    )
end

function build_expected_richness_table(df_system::DataFrame, df_community::DataFrame)
    rows = NamedTuple[]

    for row in eachrow(df_system)
        curve = compute_expected_richness_curve(
            df_system,
            df_community;
            S=row.S,
            seed=row.seed,
        )
        push!(rows, curve)
    end

    return DataFrame(rows)
end

function plot_expected_richness_by_seed(
    df_system::DataFrame,
    df_community::DataFrame;
    S::Int,
    seeds::Union{Nothing, Vector{Int}}=nothing,
    xmin::Float64=1.0,
    ymax=nothing,
    linewidth::Float64=1.2,
    linealpha::Float64=0.35,
    color::Any="#7C9FB0",
    plot_size::Tuple{Int, Int}=(520, 360),
)
    curves = build_expected_richness_table(df_system, df_community)
    curves = subset(curves, :S => ByRow(isequal(S)))
    if isnothing(seeds)
        seeds = sort(unique(curves.seed))
    end
    curves = subset(curves, :seed => ByRow(in(Set(seeds))))

    xmax = maximum(last(q) for q in curves.Qs)
    ymax_val = isnothing(ymax) ? maximum(maximum(y) for y in curves.y) : ymax

    plt = plot(
        size=plot_size,
        legend=false,
        grid=false,
        xaxis=:log10,
        xlims=(xmin, xmax),
        ylims=(0.0, ymax_val),
        xlabel="Q",
        ylabel="E[Richness | Q]",
        title="Average richness vs. Q, S=$(S)",
        framestyle=:box,
    )

    for row in eachrow(curves)
        plot!(
            plt,
            row.Qs,
            row.y;
            color=color,
            alpha=linealpha,
            linewidth=linewidth,
            label=nothing,
        )
    end

    return plt
end

function plot_expected_richness_mean_std(
    df_system::DataFrame,
    df_community::DataFrame;
    S::Int,
    seeds::Union{Nothing, Vector{Int}}=nothing,
    xmin::Float64=1.0,
    ymax=nothing,
    xmax::Float64=150.0,
    n_points::Int=300,
    seed_linewidth::Float64=0.95,
    mean_linewidth::Float64=2.8,
    mean_underlay_linewidth::Float64=5.2,
    seed_alpha::Float64=0.20,
    band_alpha::Float64=0.12,
    seed_color::Any="#92A9BD",
    mean_underlay_color::Any="#E9F1F7",
    mean_color::Any="#285A84",
    plot_size::Tuple{Int, Int}=(520, 360),
    show_seeds::Bool=true,
)
    curves = build_expected_richness_table(df_system, df_community)
    curves = subset(curves, :S => ByRow(isequal(S)))
    if isnothing(seeds)
        seeds = sort(unique(curves.seed))
    end
    curves = subset(curves, :seed => ByRow(in(Set(seeds))))

    rows = collect(eachrow(curves))
    stats = estimate_curve_mean_std(rows; x_range=(xmin, xmax), n_points=n_points)
    isnothing(stats) && error("No valid richness curves found for S=$(S)")

    ymax_data = show_seeds ? maximum(maximum(y) for y in curves.y) : maximum(stats.upper)
    ymax_val = isnothing(ymax) ? ymax_data : ymax

    plt = plot(
        size=plot_size,
        legend=false,
        grid=false,
        xaxis=:log10,
        xlims=(xmin, xmax),
        ylims=(0.0, ymax_val),
        xlabel="Q",
        ylabel="Richness",
        title="Average richness, S=$(S)",
        framestyle=:box,
    )

    if show_seeds
        for row in rows
            plot!(
                plt,
                row.Qs,
                row.y;
                color=seed_color,
                alpha=seed_alpha,
                linewidth=seed_linewidth,
                label=nothing,
            )
        end
    end

    plot!(
        plt,
        stats.xs,
        stats.lower;
        fillrange=stats.upper,
        fillalpha=band_alpha,
        linewidth=0,
        color=mean_color,
        label=nothing,
    )
    plot!(
        plt,
        stats.xs,
        stats.mean;
        color=mean_underlay_color,
        linewidth=mean_underlay_linewidth,
        label=nothing,
    )
    plot!(
        plt,
        stats.xs,
        stats.mean;
        color=mean_color,
        linewidth=mean_linewidth,
        label=nothing,
    )

    return plt
end

function build_fig4_data(;
    dir::String="data/revision",
    S::Int=6,
    xmin::Float64=1.0,
    xmax::Float64=150.0,
    pm_ymax::Float64=0.125,
    richness_ymax=nothing,
)
    results = collect_all_results(dir)
    pm_curves = build_pm_curve_table(results.df_system, results.df_community)
    pm_curves = subset(pm_curves, :S => ByRow(isequal(S)))
    richness_curves = build_expected_richness_table(results.df_system, results.df_community)
    richness_curves = subset(richness_curves, :S => ByRow(isequal(S)))

    sizes = sort(unique(pm_curves.comm_size))
    size_colors = community_colors(string.(sizes))

    pm_rows = Any[]
    pm_groups = Any[]
    for comm_size in sizes
        sub = subset(pm_curves, :comm_size => ByRow(isequal(comm_size)))
        n_types = length(unique(sub.comm_key))
        push!(pm_groups, Dict(
            "comm_size" => comm_size,
            "n_types" => n_types,
            "color" => size_colors[string(comm_size)],
        ))
        for row in eachrow(sub)
            push!(pm_rows, Dict(
                "seed" => row.seed,
                "S" => row.S,
                "comm_size" => row.comm_size,
                "comm_key" => row.comm_key,
                "Qs" => row.Qs,
                "y" => row.y,
            ))
        end
    end

    richness_rows = Any[
        Dict(
            "seed" => row.seed,
            "S" => row.S,
            "Qs" => row.Qs,
            "y" => row.y,
        )
        for row in eachrow(richness_curves)
    ]

    if isnothing(richness_ymax)
        richness_ymax = maximum(maximum(row["y"]) for row in richness_rows)
    end

    return Dict(
        "meta" => Dict(
            "S" => S,
            "x_range" => [xmin, xmax],
            "pm_ymax" => pm_ymax,
            "richness_ymax" => richness_ymax,
        ),
        "pm_groups" => pm_groups,
        "pm_rows" => pm_rows,
        "richness_rows" => richness_rows,
    )
end

function export_fig4_data(
    outpath::String="data/output/fig4_data.json";
    dir::String="data/revision",
    S::Int=6,
    xmin::Float64=1.0,
    xmax::Float64=150.0,
    pm_ymax::Float64=0.125,
    richness_ymax=nothing,
)
    data = build_fig4_data(
        ;
        dir=dir,
        S=S,
        xmin=xmin,
        xmax=xmax,
        pm_ymax=pm_ymax,
        richness_ymax=richness_ymax,
    )
    mkpath(dirname(outpath))
    open(outpath, "w") do io
        JSON.print(io, data)
    end
    return outpath
end

function export_fig4_si_data(;
    outdir::String="data/output",
    dir::String="data/revision",
    S_list=[2, 4, 8],
    xmin::Float64=1.0,
    xmax::Float64=150.0,
    pm_ymax::Float64=0.125,
    richness_ymax=nothing,
)
    mkpath(outdir)
    outpaths = String[]
    for S in S_list
        outpath = joinpath(outdir, "fig4_data-S$(S).json")
        export_fig4_data(
            outpath;
            dir=dir,
            S=S,
            xmin=xmin,
            xmax=xmax,
            pm_ymax=pm_ymax,
            richness_ymax=richness_ymax,
        )
        push!(outpaths, outpath)
    end
    return outpaths
end

if abspath(PROGRAM_FILE) == @__FILE__
    export_fig4_si_data()
end
