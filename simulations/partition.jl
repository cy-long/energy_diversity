""" Partition-specific functions: enumerate communities, compute feasibility volumes, and visualize partitions. """

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

function community_label_compact(comm)
    active = findall(comm.mask)
    isempty(active) && return "∅"
    maximum(active) <= 9 ? join(active) : "{" * join(active, ",") * "}"
end

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
        colors=colors, lws=lws,
    )
end

function build_seed_partition_run(
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
        data = build_partition_geometry_family(
            p;
            Q_range=select_range(p),
            include_empty=include_empty,
            sampling_chains=sampling_chains,
            n_sample=n_sample,
            n_layer=n_layer,
            save_centers=save_centers,
            show_prog=show_prog,
            prog_dt=prog_dt,
        )
        systems[i] = (
            seed=seed,
            S=p.S,
            σsc=σsc,
            d0=d0,
            N0=N0,
            Q_range=data.Q_range,
            vols_C=data.vols_C,
            communities=copy(data.community_table),
            centers=save_centers ? data.all_centers : nothing,
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
