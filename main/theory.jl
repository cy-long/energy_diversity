""" Analyze the results in theory """

using Statistics, StatsPlots, StatsBase, Combinatorics, CategoricalArrays, Distributions
using DataFrames, Glob, JLD2
using Plots, Interpolations

# read all results at data/main
function collect_all_results(dir::String)
    files = glob("results_seed*.jld2", dir)
    all_rows = DataFrame[]
    
    for file in files
        @info "Loading $file"
        data = load(file, "results")  # Vector{Dict}

        # Extract seed from filename
        seed_match = match(r"seed(\d+)", file)
        seed = isnothing(seed_match) ? missing : parse(Int, seed_match.captures[1])

        # Add :seed key to each dictionary
        rows = [merge(d, Dict(:seed => seed)) for d in data]
        push!(all_rows, DataFrame(rows))
    end
    
    col_orders = [:seed, :S, :σsc, :d0, :N0, :type, :Qs, :vols, :devols]
    combined = vcat(all_rows...)
    return combined[!, col_orders]
end

function estimate_function_mean(df::DataFrame; x_range::Tuple, n_points::Int=300)
    xs = 10 .^ range(log10(x_range[1]), log10(x_range[2]), length=n_points)
    y_interp_stack = Matrix{Float64}(undef, length(xs), 0)

    for row in eachrow(df)
        x = row.Qs
        y = row.vols ./ row.devols
        itp = LinearInterpolation(x, y, extrapolation_bc=Flat()) # Linear, flat extrap
        y_interp = itp.(xs)
        y_interp_stack = hcat(y_interp_stack, y_interp)
    end

    ϵ = 1e-4  # small positive floor
    log_y = log.(clamp.(y_interp_stack, ϵ, 1.0))
    log_mean = vec(mean(log_y; dims=2))
    log_std = vec(std(log_y; dims=2))

    return xs, log_mean, log_std
end

function show_single(df::DataFrame, Q_range::Tuple=(1.0,100.0); p_size::Tuple=(160,120), type::String, c_params::Tuple=(1.0, 1.0, 1.0), y_lim::Union{Float64, Nothing}=nothing)

    legends = Dict("r"=>(0.11,0.93), "i"=> (0.11,0.32));
    p_kwargs = Dict(
        :size => p_size,
        :framestyle => :box,
        :xlabel => "Total Energy Supply",
        :ylabel => labels[type],
        :grid => false,
        :xaxis => :log10,
        :xlim => Q_range,
        :guidefont => font("Helvetica", 6),
        :tickfont => font("Helvetica", 5),
        :legendfont => font("Helvetica", 5),
        :legend=> legends[type],
        :foreground_color_legend => nothing,
        :legend_border => false,
        :background_color => :transparent,
        :linewidth => 1.0,
        :yticks => [0.25,0.75]
    )
    if y_lim !== nothing 
        p_kwargs[:ylim] = (0.0, y_lim)
    end
    plt = plot(; p_kwargs...)

    S = 2; cl = colors[type][1]
    df_S = filter(row -> row.S==S && row.σsc==c_params[1] && row.d0==c_params[2] && row.N0==c_params[3], df);

    if df_S == DataFrame()
        @warn "No data for S=$(S), σsc=$(c_params[1]), d0=$(c_params[2]), N0=$(c_params[3])"
        return nothing
    end

    row =  eachrow(df_S)[25]
        plot!(plt, row.Qs, row.vols ./ row.devols, alpha=0.5, lw=1.0, color = cl, label="")
    return plt
end

function show_group_curves(df::DataFrame, Q_range::Tuple=(1.0,100.0); p_size::Tuple=(240,200), type::String, c_params::Tuple=(1.0, 1.0, 1.0), y_lim::Union{Float64, Nothing}=nothing)
    legends = Dict("r"=>(0.11,0.93), "i"=> (0.11,0.32));
    p_kwargs = Dict(
        :size => p_size,
        :framestyle => :box,
        :xlabel => "Total Energy Supply",
        :ylabel => labels[type],
        :grid => false,
        :xaxis => :log10,
        :xlim => Q_range,
        :guidefont => font("Helvetica", 7),
        :tickfont => font("Helvetica", 6),
        :legendfont => font("Helvetica", 6),
        :legend=> legends[type],
        :foreground_color_legend => nothing,
        :legend_border => false,
        :background_color => :transparent,
        :linewidth => 1.0,
    )
    if y_lim !== nothing 
        p_kwargs[:ylim] = (0.0, y_lim)
    end
    plt = plot(; p_kwargs...)

    for (S, cl) in zip([8,6,4,2], reverse(colors[type]))
        df_S = filter(row -> row.S==S && row.σsc==c_params[1] && row.d0==c_params[2] && row.N0==c_params[3], df);
        if df_S == DataFrame()
            @warn "No data for S=$(S), σsc=$(c_params[1]), d0=$(c_params[2]), N0=$(c_params[3])"
            return nothing
        end
        xs, log_mean, log_std = estimate_function_mean(df_S, x_range=Q_range);
        lower = exp.(log_mean .- log_std); upper = exp.(log_mean .+ log_std);
        plot!(plt, xs, lower, fillrange=upper, fillalpha=0.35, color=cl, label="", linealpha=0)
        plot!(plt, xs, exp.(log_mean), color=cl, xaxis=:log10, label="")
    end
    for (S, cl) in zip([2,4,6,8], colors[type])
        plot!(plt, [NaN], [NaN], color=cl, label="S=$(S)")
    end
    return plt
end

function show_inset(df::DataFrame; type::String, c_params::Tuple=(1.0, 1.0, 1.0))
    colors = Dict("r"=> :midnightblue, "i"=> :firebrick);
    yticks = Dict("r" => [0.01], "i" => [0.025]);
    ylims = Dict("r" => (0.0, 0.015), "i" => (0.0, 0.05));
    df_S8 = filter(row -> row.S==8 && row.σsc==c_params[1] && row.d0==c_params[2] && row.N0==c_params[3], df);
    
    if df_S8 != DataFrame()
        xs, log_mean, log_std = estimate_function_mean(df_S8, x_range=(1,100));
        lower = exp.(log_mean .- log_std); upper = exp.(log_mean .+ log_std);
        insetplt = plot(
            xs, exp.(log_mean), 
            color=colors[type], 
            linewidth=0.8, 
            xaxis=:log10, 
            legend=false, 
            size=(90,60), 
            xlim=(8,100), 
            xlabel="", 
            ylabel="", 
            xticks=[10,100], 
            ylim=ylims[type], 
            yticks=yticks[type],
            tickfont = font("Helvetica", 3),
            grid = false,
            foreground_color_legend = nothing,
            background_color = :transparent,
        )
        plot!(insetplt, xs, lower, fillrange=upper, fillalpha=0.35, color=colors[type], linealpha=0)
        return insetplt
    end
end


results_df = collect_all_results("data/main")
# @load "data/results_df.jld2" results_df

results_r = results_df[results_df.type .== :total, :];
results_i = results_df[results_df.type .== :indiv, :];

colors = Dict(
    "r"=> [:lightblue, :dodgerblue, :royalblue, :midnightblue], 
    "i"=> [:lightsalmon, :darkorange, :orangered, :firebrick]
);

labels = Dict(
    "r" => "Prob. Maturation",
    "i" => "Prob. Developing",
);

# --- general patterns --- S=4, all scale=1; unimodal vs. saturation
## caveat: 1 std tends to underestimate the variance, is just a way to visualize the uncertainty

plt1 = show_group_curves(results_r, type="r")
# savefig(plt1, "figures/curves_r.svg")
plt1i = show_inset(results_r, type = "r")
# savefig(plt1i, "figures/inset_r.svg")

plt2 = show_group_curves(results_i, type="i")
# savefig(plt2, "figures/curves_i.svg")
plt2i = show_inset(results_i, type = "i")
# savefig(plt2i, "figures/inset_i.svg")

plt3 = show_single(results_r, type="r",y_lim=1.0)
# savefig(plt3, "figures/methods_r2.svg")

# --- scores ---
df = results_r; c_params = (1.0, 1.0, 1.0); 
Gss = Vector{Vector{Float64}}();
for S in [2,4,6,8]
    i = Int(S/2)
    df_S = filter(row -> row.S==S && row.σsc==c_params[1] && row.d0==c_params[2] && row.N0==c_params[3], df);
    Gs = Vector{Float64}();
    for row in eachrow(df_S)
        push!(Gs, score_unimodal(row.Qs, row.vols ./ row.devols))
    end
    push!(Gss, Gs)
end
mean(Gss[1])
mean(Gss[2])
mean(Gss[3])
mean(Gss[4])