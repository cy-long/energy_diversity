using EnerFeas
using Plots
using DataFrames

function create_color(names::Vector{String})
    colors = Dict{String, Any}()
    palette = [
        :steelblue3, :indianred2, :darkorange2, :seagreen3, :mediumpurple3,
        :goldenrod2, :teal, :firebrick3, :dodgerblue3, :sienna3
    ]
    L = length(palette)
    for (i, name) in enumerate(names)
        colors[name] = palette[mod1(i - 1, L)]
    end
    colors[names[1]] = :gray80
    return colors
end

function create_alpha(names::Vector{String})
    N = length(names)
    S = log2(N)
    alphas = Dict{String, Float64}()
    for name in names
        k = length(name)
        alphas[name] = 0.5 + 0.5 * (k / S)
    end
    return alphas
end

function create_lw(names::Vector{String})
    N = length(names)
    S = log2(N)
    lws = Dict{String, Float64}()
    for name in names
        k = length(name)
        lws[name] = 0.5 + 2 * (k / S)
    end
    return lws
end

function _community_label_compact(comm)
    active = findall(comm.mask)
    isempty(active) && return "∅"
    if maximum(active) <= 9
        return join(active)
    end
    return "{" * join(active, ",") * "}"
end

community_label_compact(comm) = _community_label_compact(comm)

function enumerate_comm(S::Int; include_empty::Bool=true)
    comms = enumerate_communities(S; include_empty=include_empty)
    names = [community_label_compact(comm) for comm in comms]
    return comms, names
end

function build_partition_geometry_family(
    p;
    Q_range::Vector{Float64}=select_range(p),
    include_empty::Bool=true,
    sampling_chains::Int=2,
    n_sample::Int=5 * 10^4,
    n_layer::Int=10,
    save_centers::Bool=false,
    show_prog::Bool=true,
    prog_dt::Float64=0.5,
)
    comms, names = enumerate_comm(p.S; include_empty=include_empty)
    N = length(comms)

    all_volumes = Dict{String, Vector{Float64}}()
    all_centers = save_centers ? Dict{String, Vector}() : nothing
    colors = create_color(names)
    alphas = create_alpha(names)
    lws = create_lw(names)
    community_rows = NamedTuple[]

    vols_C = volume_range_C(p, Q_range)
    for (name, comm) in zip(names, comms)
        @info "Probe EFD of community $name"
        vols_c = if save_centers
            vols_tmp, cens_c = compute_range_EFD(
                p, Q_range, comm;
                n_chains=sampling_chains, n_sample=n_sample, n_layer=n_layer,
                show_prog=show_prog, prog_dt=prog_dt
            )
            all_centers[name] = cens_c
            vols_tmp
        else
            volume_range_EFD(
                p, :matr, Q_range, comm;
                n_chains=sampling_chains, n_sample=n_sample, n_layer=n_layer,
                show_prog=show_prog, prog_dt=prog_dt
            )
        end
        all_volumes[name] = vols_c
        push!(community_rows, (
            comm_active=collect(comm.active),
            comm_key=name,
            comm_size=length(comm.active),
            volumes=vols_c,
            pm=vols_c ./ vols_C,
        ))
    end

    return (
        p=p, S=p.S, N=N, Q_range=Q_range,
        comms=comms, names=names,
        all_volumes=all_volumes, all_centers=all_centers, vols_C=vols_C,
        community_table=DataFrame(community_rows),
        colors=colors, alphas=alphas, lws=lws,
    )
end

function build_system_partition_block(
    p;
    seed::Int,
    σsc::Float64=1.0,
    d0::Float64=1.0,
    N0::Float64=1.0,
    Q_range::Vector{Float64}=select_range(p),
    include_empty::Bool=true,
    sampling_chains::Int=2,
    n_sample::Int=5 * 10^4,
    n_layer::Int=10,
    save_centers::Bool=false,
    show_prog::Bool=true,
    prog_dt::Float64=0.5,
)
    data = build_partition_geometry_family(
        p;
        Q_range=Q_range,
        include_empty=include_empty,
        sampling_chains=sampling_chains,
        n_sample=n_sample,
        n_layer=n_layer,
        save_centers=save_centers,
        show_prog=show_prog,
        prog_dt=prog_dt,
    )

    communities = copy(data.community_table)
    return (
        seed=seed,
        S=p.S,
        σsc=σsc,
        d0=d0,
        N0=N0,
        Q_range=data.Q_range,
        vols_C=data.vols_C,
        communities=communities,
        centers=save_centers ? data.all_centers : nothing,
    )
end

function build_seed_partition_survey(
    seed::Int;
    S_range::Vector{Int}=[2, 4, 6, 8],
    σsc::Float64=1.0,
    d0::Float64=1.0,
    N0::Float64=1.0,
    include_empty::Bool=true,
    sampling_chains::Int=2,
    n_sample::Int=5 * 10^4,
    n_layer::Int=10,
    save_centers::Bool=false,
    show_prog::Bool=true,
    prog_dt::Float64=0.5,
)
    grand_S = maximum(S_range)
    p0 = generate_model_system(grand_S, seed, σsc, d0, N0)
    systems = Vector{NamedTuple}(undef, length(S_range))

    for (i, S) in pairs(S_range)
        p = sub_model_system(S, p0)
        systems[i] = build_system_partition_block(
            p;
            seed=seed,
            σsc=σsc,
            d0=d0,
            N0=N0,
            Q_range=select_range(p),
            include_empty=include_empty,
            sampling_chains=sampling_chains,
            n_sample=n_sample,
            n_layer=n_layer,
            save_centers=save_centers,
            show_prog=show_prog,
            prog_dt=prog_dt,
        )
    end

    return (
        seed=seed,
        σsc=σsc,
        d0=d0,
        N0=N0,
        S_range=collect(S_range),
        systems=systems,
    )
end

function collect_partition_tables(results)
    df_system = DataFrame(
        seed=Int[],
        S=Int[],
        σsc=Float64[],
        d0=Float64[],
        N0=Float64[],
        Q_range=Vector{Float64}[],
        vols_C=Vector{Float64}[],
    )
    community_tables = DataFrame[]

    for block in results.systems
        push!(df_system, (
            seed=block.seed,
            S=block.S,
            σsc=block.σsc,
            d0=block.d0,
            N0=block.N0,
            Q_range=block.Q_range,
            vols_C=block.vols_C,
        ))

        table = copy(block.communities)
        insertcols!(table, 1, :S => fill(block.S, nrow(table)))
        insertcols!(table, 1, :seed => fill(block.seed, nrow(table)))
        push!(community_tables, table)
    end

    df_community = isempty(community_tables) ? DataFrame() : vcat(community_tables...)
    return df_system, df_community
end


function compute_prob_maturation_curves(data; empty_name::String="∅")
    pm_by_name = Dict{String, Vector{Float64}}()
    for name in data.names
        pm_by_name[name] = data.all_volumes[name] ./ data.vols_C
    end

    nonempty_names = [name for name in data.names if name != empty_name]
    pm_nonempty = zeros(length(data.Q_range))
    for name in nonempty_names
        pm_nonempty .+= pm_by_name[name]
    end

    return (pm_by_name=pm_by_name, pm_nonempty=pm_nonempty, nonempty_names=nonempty_names)
end

function plot_prob_maturation_curves(data; empty_name::String="∅")
    pm = compute_prob_maturation_curves(data; empty_name=empty_name)
    plt = plot(legend=:topright, grid=false, size=(480, 374), legend_foreground_color=nothing)
    plot!(plt, [], []; label=empty_name, color=:gray80, linewidth=0.1)
    for name in data.names
        name == empty_name && continue
        plot!(plt, data.Q_range, pm.pm_by_name[name];
            label=name, linewidth=data.lws[name], color=data.colors[name])
    end
    plot!(plt, data.Q_range, pm.pm_nonempty;
        label="non-$empty_name", linewidth=2, color=:black)
    xlabel!(plt, "Q")
    ylabel!(plt, "Pᴹ")
    xlims!(plt, (data.Q_range[1], data.Q_range[end]))
    ylims!(plt, (0.0, 1.0))
    xaxis!(plt, :log10)
    return plt
end

function plot_partition_volumes(data)
    plt = plot(legend=:outerright, grid=false, size=(480, 374), legend_foreground_color=nothing)
    for name in data.names
        plot!(plt, data.Q_range, data.all_volumes[name];
            linewidth=data.lws[name], color=data.colors[name], label=name)
    end
    xlabel!(plt, "Q")
    ylabel!(plt, "vol")
    xlims!(plt, (data.Q_range[1], data.Q_range[end]))
    ylims!(plt, (1e-3, 1e3))
    xaxis!(plt, :log10)
    yaxis!(plt, :log10)
    yticks!(plt, [1e-3, 1e0, 1e3])
    return plt
end

function compute_average_richness_curve(data)
    Nq = length(data.Q_range)
    R = zeros(Float64, Nq)

    for (name, comm) in zip(data.names, data.comms)
        r_i = count(comm.mask)
        P_i = data.all_volumes[name] ./ data.vols_C
        R .+= r_i .* P_i
    end

    return (Q_range=data.Q_range, R=R)
end

function plot_average_richness_curve(data)
    rc = compute_average_richness_curve(data)
    plt = plot(legend=false, grid=false, size=(480, 374))
    plot!(plt, rc.Q_range, rc.R; linewidth=2.5, color=:black)
    xlabel!(plt, "Q")
    ylabel!(plt, "R(Q)")
    xlims!(plt, (rc.Q_range[1], rc.Q_range[end]))
    xaxis!(plt, :log10)
    return plt
end
