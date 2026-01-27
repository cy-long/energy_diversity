""" A case study on human gut microbiome and host nutrient intake"""

using CSV, DataFrames, Plots, Random
using Statistics, StatsPlots, StatsBase, Combinatorics, CategoricalArrays
using ProgressMeter
using Printf
using ColorSchemes

function filter_taxa(df::DataFrame, thre::Float64)
    remove_taxa = String[];
    for col in names(df)[3:end]
        vals = df[!, col]
        if col == "" || all(<=(thre), vals) || all(>=(thre), vals)
            push!(remove_taxa, col)
        end
    end
    final_df = select(df, Not(remove_taxa));
    final_taxa = names(final_df)[3:end];
    @info "$(length(final_taxa)) final candidates under threshold $thre \n"
    return final_df, final_taxa
end

function categorize_calories(df::DataFrame)
    calories = df.calories;
    cals_cut = quantile(calories[.!isnan.(calories)], [1/3, 2/3]);
    @info "Calorie cutoffs: $(cals_cut[1]), $(cals_cut[2]) in $(length(calories[.!isnan.(calories)])) days)"
    # if cals_cut === nothing
    #     n_cal = length(df.calories[.!isnan.(df.calories)])
    #     @info "Determining calorie cutoffs using $n_cal days of info)"
    #     cals_cut = quantile(df.calories[.!isnan.(df.calories)], [1/3, 2/3]);
    # end
    levels = ["low", "medium", "high"]
    df.calories = ifelse.(df.calories .< cals_cut[1], levels[1],
        ifelse.(df.calories .<= cals_cut[2], levels[2], levels[3]));
    df.calories = categorical(df.calories, levels=levels, ordered=true);
    counts = countmap(df.calories);
    return df, counts
end

function summarize_community(df::DataFrame, community::Vector{Symbol}, cal_counts::Dict, thre::Float64)
    df_c = select(df, [:calories, community...])
    f_levels = Float64[]; occurrence_ = 0; levels = ["low", "medium", "high"]
    cal_count_unified = maximum(values(cal_counts))
    for level in levels
        occur_count = sum(all(row .> thre) for row in eachrow(Matrix(df_c[df_c.calories .== level, community])))
        push!(f_levels, occur_count /cal_count_unified) # avoid artifact from different calorie counts!! since this is 1/3 2/3 cut
        occurrence_ += occur_count
    end
    if occurrence_ == 0
        return nothing # skip communities that never co-present
    end
    comm_name = join([string(sym)[1:3] for sym in community], ",")
    idx = findall(x -> x == maximum(f_levels), f_levels)
    # _, idx = findmax(f_levels)
    # if length(idx) > 1
    #     peak_idx = "tie"
    # else
    #     peak_idx = levels[idx[1]]
    # end
    return (community=comm_name, freq_low=f_levels[1], freq_medium=f_levels[2], freq_high=f_levels[3], occurrence=occurrence_, peak=idx);
end


function summarize_size(df::DataFrame, taxa::Vector{String}, K::Int, thre::Float64, cal_counts::Dict)
    communities = collect(combinations(Symbol.(taxa), K));
    results = DataFrame(
        community = String[],
        freq_low = Float64[], 
        freq_medium = Float64[], 
        freq_high = Float64[], 
        occurrence = Int[], 
        peak = Vector{Int}[]
    )
    for comm in communities
        r = summarize_community(df, comm, cal_counts, thre)
        r !== nothing && push!(results, r)
    end
    return results
end

function summarize_system(df, taxa, K_ranges, thre, cal_counts, occ_cut)
    all_counts = DataFrame();
    level_map = Dict(1 => "low", 2 => "medium", 3 => "high")
    cut_labels = ["<=$occ_cut", ">$occ_cut"]
    @showprogress for K in K_ranges
        results = summarize_size(df, taxa, K, thre, cal_counts);
        results.K = fill(K, nrow(results))
        results.occur_bin = ifelse.(results.occurrence .<=occ_cut, cut_labels[1], cut_labels[2])
        results.occur_bin = categorical(results.occur_bin; levels=cut_labels, ordered=true)

        # Flatten peaks with proportional weights
        rows = []
        for row in eachrow(results)
            if isempty(row.peak)
                push!(rows, (K = row.K, occur_bin = row.occur_bin, peak = "tie", weight = 1.0))
            else
                w = 1.0 / length(row.peak)
                for pk in row.peak
                    peak_str = haskey(level_map, pk) ? level_map[pk] : "tie"
                    push!(rows, (K = row.K, occur_bin = row.occur_bin, peak = peak_str, weight = w))
                end
            end
        end
        flat_df = DataFrame(rows)

        counts_K = combine(groupby(flat_df, [:K, :occur_bin, :peak]), :weight => sum => :count)
        counts_K = combine(groupby(counts_K, [:K, :occur_bin])) do df
            total = sum(df.count)
            df.ratio = df.count ./ total
            df
        end

        counts_K.K = fill(string(K), nrow(counts_K)) # Store K as string for x-axis grouping
        append!(all_counts, counts_K)
    end
    all_counts.K = categorical(all_counts.K, ordered=true);
    all_counts.occur_bin = categorical(all_counts.occur_bin, levels=["<=1", ">1"], ordered=true);
    all_counts.peak = categorical(all_counts.peak, levels=["high", "medium", "low", "tie"], ordered=true);
    return all_counts
end

function show_rep_communities(results::DataFrame, top::Int, ymax=0.25)
    top_comms = sort(results, :occurrence, rev=true)[1:top,:];
    top_comms_long = stack(
        top_comms,
        [:freq_low, :freq_medium, :freq_high],
        variable_name = :calorie,
        value_name = :frequency
    )
    top_comms_long[!, :type] .= "m"

    bottom_comms = sort(results, :occurrence, rev=false)[1:top,:];
    bottom_comms_long = stack(
        bottom_comms,
        [:freq_low, :freq_medium, :freq_high],
        variable_name = :calorie,
        value_name = :frequency
    )
    bottom_comms_long[!, :type] .= "d"

    rep_comms_long = vcat(top_comms_long, bottom_comms_long);
    rep_comms_long.calorie = replace.(rep_comms_long.calorie, r"freq_" => "")

    colors = Dict("m"=>:dodgerblue, "d"=>:darkorange);
    linestyles = Dict("m"=>:solid, "d"=>:dash);
    markershapes = Dict("m"=>:circle, "d"=>:rect);

    pl = plot(size=(160,100), 
        grid=false, legend=:top, linewidth = 1.0,
        guidefont = font("Helvetica", 6), tickfont=font("Helvetica", 5),
        labelfont = font("Helvetica", 5), titlefont=font("Helvetica", 7),
        foreground_color_legend=nothing,
        background_color=:transparent,
        legend_border=false,
        markercolor=:white,
    )

    @df rep_comms_long plot!(pl,
        :calorie,
        :frequency,
        group = :community,
        xlabel = "",
        ylabel = "",
        color = [colors[t] for t in :type],
        label = "",
        ylim=(-0.06, ymax),
        yticks = [0.0, 0.4, 0.8],
        markercolor= :white,
        linestyle = [linestyles[t] for t in :type],
        markershape = [markershapes[t] for t in :type],
        markerstrokecolor = [colors[t] for t in :type],
        markersize=2,
        markerstrokewidth = 1,
    )

    # overlay peak frequencies
    calorie_idx = Dict("low"=>1, "medium"=>2, "high"=>3);

    peak_rows = filter(row -> calorie_idx[row.calorie] in row.peak, eachrow(rep_comms_long))
    peak_df = DataFrame(peak_rows)

    @df peak_df scatter!(
        pl,
        :calorie,
        :frequency,
        label = "",
        color = [colors[t] for t in :type],
        markershape = [markershapes[t] for t in :type],
        markercolor = [colors[t] for t in :type],
        markerstrokecolor = [colors[t] for t in :type],
        markersize = 2,
        markerstrokewidth = 2,
    )

    # pll = plot(
    #     size=(120,90), 
    #     grid=false, legend=:top,
    #     guidefont = font("Helvetica", 6), tickfont=font("Helvetica", 5), labelfont = font("Helvetica", 5), titlefont=font("Helvetica", 6),
    #     foreground_color_legend=nothing,
    #     background_color=:transparent,
    #     legend_border=false,
    #     markersize=2,
    #     markercolor=:white,
    #     framestyle=:none,
    #     xaxis=false, yaxis=false,
    #     xticks=false, yticks=false,
    # )

    # plot!(pll, [NaN], [NaN], label="most occurring (maturation)", color=:dodgerblue, linestyle=:solid, markershape = :circle, markercolor=:white, markerstrokecolor = :dodgerblue, markersize=5, markerstrokewidth = 2)
    # plot!(pll, [NaN], [NaN], label="least occurring (developing)", color=:darkorange, linestyle=:solid, markershape = :rect, 
    # markercolor=:white, markerstrokecolor = :darkorange, markersize=5, markerstrokewidth = 2)
    # scatter!(pll, [NaN], [NaN], label="peak", color=:dodgerblue, markershape=:circle, markerstrokecolor = :dodgerblue, markersize=4, markerstrokewidth = 2,)
    # scatter!(pll, [NaN], [NaN], label="peak", color=:darkorange, markershape=:rect, markerstrokecolor = :darkorange, markersize=4, markerstrokewidth = 2,)

    # return pl, pll
    return pl
end

function show_all_counts(all_counts, size=(160,100))
    unique_K = sort!(unique(all_counts.K));
    occur_bins = ["<=1", ">1"];
    xvals = Dict{Tuple{String, String}, Float64}();
    offset = 0.0
    spacing_within = 0.4
    spacing_between = 1.0
    x_ticks = Float64[]
    x_tick_labels = String[]

    for k in unique_K
        for (i, bin) in enumerate(occur_bins)
            key = (string(k), bin)
            x = offset + (i-1) * spacing_within
            xvals[key] = x
            push!(x_ticks, x)
            push!(x_tick_labels, "$(k)-$bin")
        end
        offset += spacing_between
    end

    all_counts.x = [xvals[(string(row.K), row.occur_bin)] for row in eachrow(all_counts)];
    xtick_positions = [mean([xvals[(string(k), b)] for b in occur_bins]) for k in unique_K];
    xtick_labels = String.(unique_K);

    Cmap = Dict(
        "low" => :khaki4,
        "medium" => :rosybrown,
        "high" => :darkolivegreen
    )
    bar_colors = [Cmap[string(row.peak)] for row in eachrow(all_counts)];
    fillalphas = ifelse.(all_counts.occur_bin .== "<=1", 0.1, 1.0)

    plt = plot(size=size, grid=false);
    @df all_counts groupedbar!(
        plt,
        :x,
        :ratio,
        group = :peak,
        bar_position = :stack,
        xlabel = "",
        ylabel = "",
        guidefont = font("Helvetica", 6), tickfont=font("Helvetica", 5),
        labelfont = font("Helvetica", 5), titlefont=font("Helvetica", 7),
        foreground_color_legend=nothing,
        background_color=:transparent,
        bar_width = 0.3,
        label="",
        xticks = (xtick_positions, xtick_labels),
        linealpha = 1.0,
        linewidth = 1.0,
        color = bar_colors,
        linecolor = bar_colors,
        fillalpha = fillalphas,
        yticks=[0.0, 0.5, 1.0],
        ytick_direction=:out,
    );

    for peak in ["high", "medium", "low"] # change orders to high mid low
        col = Cmap[peak]
        bar!(plt, [NaN], [NaN], color=col, alpha=1.0, label="$(peak)", bar_width=0.3, linecolor = col)
    end

    bar!(plt, [NaN], [NaN], color=:white, alpha=0.0, label="≤1", bar_width=0.3, linecolor = :gray10)
    bar!(plt, [NaN], [NaN], color=:gray10, alpha=1.0, label=">1", bar_width=0.3, linecolor = :gray10)

    plot!(plt, legend=:outerright, legendfont=("Helvetica", 4),
        foreground_color_legend=nothing, legend_border=false
    );

    return(plt)
end

# -------------- Processing the dataset --------------
# --- Read metadata (calories, timestamps) ---
meta_keys = CSV.read("data/lifestyle/metadata/meta.keys", DataFrame; header=false)[!,2];
metadata = CSV.read("data/lifestyle/metadata/subjectA.gut.M.proc", DataFrame; header=false);
rename!(metadata, Symbol.(meta_keys));

calories = metadata[!, :nutrition_calorie];
# fat = metadata[!, :nutrition_fat];
# protein = metadata[!, :nutrition_protein];
# carbon = metadata[!, :nutrition_carb];


timestamps = CSV.read("data/lifestyle/metadata/meta.time", DataFrame; header=false)[!,1];
time_calo = timestamps[findall(calories .> 0)];

# pl_A = plot(timestamps, calories, marker=:circle, xlim=(0, 300), markerstrokecolor=:auto, grid=false, framestyle=:box, label="", size=(180,90), markersize=1, guidefont = font("Helvetica", 6), tickfont=font("Helvetica", 5), labelfont = font("Helvetica", 5), foreground_color_legend=nothing, legend_border=false, background_color=:transparent, lw=0.5);
# pl_A
# savefig(pl_A, "figures/calorie_time.svg")

# pl_B = plot(size=(180,90), grid=false, guidefont = font("Helvetica", 6), tickfont=font("Helvetica", 5), labelfont = font("Helvetica", 5), foreground_color_legend=nothing, legend_border=false, background_color=:transparent, yticks=[4,8,12]);
# pl_B = histogram!(pl_B, calories, bins=50, title="", xlabel="", ylabel="", legend=false, grid=false, color=:lightblue, alpha=0.5);
cals_cut_plot = quantile(calories[.!isnan.(calories)], [1/3, 2/3]);
# vline!(pl_B, [cals_cut_plot[1]], color=:gray90, linestyle=:dash, label="");
# vline!(pl_B, [cals_cut_plot[2]], color=:gray60, linestyle=:dash, label="");
# display(pl_B)

# savefig(pl_B, "figures/calorie_cut.svg")

# --- Read taxanomy (from Order to Species) ---
microb_keys = CSV.read("data/lifestyle/microbiota/subjectA.gut.keys", DataFrame; header=true);
taxo = split.(microb_keys.Taxo, ';');
microb_keys.Order   = replace.(getindex.(taxo, 4), r"^o__" => "");
microb_keys.Family  = replace.(getindex.(taxo, 5), r"^f__" => "");
microb_keys.Genus   = replace.(getindex.(taxo, 6), r"^g__" => "");
microb_keys.Species = replace.(getindex.(taxo, 7), r"^s__" => "");
select!(microb_keys, Not(:Taxo, :Source));

# --- Read sequencing time series ---
microb_counts = CSV.read("data/lifestyle/microbiota/subjectA.gut.micro", DataFrame; header=false);

# rename columns (time and OTUs)
rename!(microb_counts, vcat(:time, Symbol.(microb_keys.OTU)));
microb_counts.time = Int.(microb_counts.time);

# filter rows (days) with calories records and combine them
# time_calo = timestamps;
microb_calo = subset!(microb_counts, :time => ByRow(in(time_calo)));
microb_calo.calories = [calories[findfirst(==(t), timestamps)] for t in microb_calo.time];
microb_calo = select(microb_calo, :time, :calories, Not([:time, :calories])); # reorder colums


# --- Aggregate OTU counts into taxonomy level ---
taxo_level = "Family" # choose from "Order", "Family", "Genus", "Species"
taxa = unique(getproperty(microb_keys, Symbol(taxo_level)));
taxa_calo = DataFrame(time = microb_calo.time, calories = microb_calo.calories);
# sum OTU counts into corresponding taxon
for taxon in taxa
    otus = Symbol.(microb_keys.OTU[microb_keys[!, Symbol(taxo_level)] .== taxon])
    taxa_calo[!, Symbol(taxon)] = sum(eachcol(select(microb_calo, otus)))
end
# remove rows (days) with that are all NaN
filter!(row -> all(!isnan, row), taxa_calo);


# --- Convert into relative abundance ---
rela_abun_calo = copy(taxa_calo);
for i in 1:nrow(rela_abun_calo)
    total = sum(rela_abun_calo[i, 3:end])
    # if total == 0
    #     continue
    # end
    for j in 3:ncol(rela_abun_calo)
        rela_abun_calo[i, j] /= total
    end
end


# pl_C = plot(
#     rela_abun_calo[!,"time"], rela_abun_calo[!,"Desulfovibrionales"],marker=:circle, markerstrokecolor=:auto, grid=false, framestyle=:box, label="", size=(180,90), markersize=1, guidefont = font("Helvetica", 6), tickfont=font("Helvetica", 5), labelfont = font("Helvetica", 5), titlefont=font("Helvetica", 6),xlim=(0,300), foreground_color_legend=nothing, legend_border=false, yticks=[0.001, 0.01],background_color=:transparent, lw=0.5, color=:gray50
# );
# hline!(pl_C, [0.001], color=:red, linestyle=:dash, label=false);
# plot!(pl_C, xlabel="", ylabel="",title="", legend=false);
# display(pl_C)
# savefig(pl_C, "figures/rela_abun_ts.svg");


# -------------- Anallyze the results --------------
# --- 1. top 5 occurring 4-taxa communities ---
K = 3
final_rela_calo, taxa_candidates = filter_taxa(rela_abun_calo, 1e-3);
final_rela_calo, cal_counts = categorize_calories(final_rela_calo);
results_K = summarize_size(final_rela_calo, taxa_candidates, K, 1e-3, cal_counts);

# pl_D = plot(size=(180,90),grid=false, label="", markersize=1, guidefont = font("Helvetica", 6), tickfont=font("Helvetica", 5), labelfont = font("Helvetica", 5), titlefont=font("Helvetica", 6), foreground_color_legend=nothing, legend_border=false,background_color=:transparent, )

# histogram!(pl_D, results_K.occurrence, normalize=false, bins=20, yaxis=:log10, ylim=(1,1000), xlim=(1,15), label="")
# plot!(pl_D, xlabel="", ylabel="",yticks=[1,100],xticks=[5,10])
# # savefig(pl_D, "figures/presence_hist_K4.svg")

pl_E = show_rep_communities(results_K, 10, 1.2)
pl_E
savefig(pl_E, "figures/most_least_occur-Fam.svg")


# ---- 2. Summarize peak distributions for K∈[2,3,4,5] ----
K_ranges = [2,3,4,5];
final_rela_calo, taxa_candidates = filter_taxa(rela_abun_calo, 1e-3);
final_rela_calo, cal_counts = categorize_calories(final_rela_calo);
taxa_candidates = sort(taxa_candidates);
all_counts = summarize_system(final_rela_calo, taxa_candidates, K_ranges, 1e-3, cal_counts, 1);

pl_F = show_all_counts(all_counts)
savefig(pl_F, "figures/peak_distributions-family.svg")