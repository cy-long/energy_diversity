using Pkg
Pkg.activate(".")
using CSV, DataFrames, Plots, Random
using Statistics, StatsPlots, StatsBase, Combinatorics, CategoricalArrays
using ProgressMeter
using Printf

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

function categorize_calories(df::DataFrame, cals_cut::Vector)
    levels = ["low", "mid", "high"]
    df.calories = ifelse.(df.calories .< cals_cut[1], levels[1],
        ifelse.(df.calories .<= cals_cut[2], levels[2], levels[3]));
    df.calories = categorical(df.calories, levels=levels, ordered=true);
    counts = countmap(df.calories);
    return df, counts
end

function summarize_community(df::DataFrame, community::Vector{Symbol}, cal_counts::Dict, thre::Float64)
    df_c = select(df, [:calories, community...])
    f_levels = Float64[]; presence = 0; levels = ["low", "mid", "high"]
    for level in levels
        presence_level = sum(all(row .> thre) for row in eachrow(Matrix(df_c[df_c.calories .== level, community])))
        push!(f_levels, presence_level / cal_counts[level])
        presence += presence_level
    end
    if presence == 0
        return nothing # skip communities that never co-present
    end
    comm_name = join([string(sym)[1:3] for sym in community], ",")
    _, idx = findmax(f_levels)
    return (community=comm_name, freq_low=f_levels[1], freq_mid=f_levels[2], freq_high=f_levels[3], presence=presence, peak=levels[idx]);
end

function show_top_communities(results::DataFrame, top::Int, K::Int, thre::Float64)
    top_comms = sort(results, :presence, rev=true)[1:top,:];
    top_comms_long = stack(
        top_comms,
        [:freq_low, :freq_mid, :freq_high],
        variable_name = :calorie,
        value_name = :frequency
    )
    top_comms_long.calorie = replace.(top_comms_long.calorie, r"freq_" => "")
    plt1 = @df top_comms_long plot(
        :calorie,
        :frequency,
        group = :community,
        linestyle = :linestyle,
        legend = :outerright,
        xlabel = "Calorie Level",
        ylabel = "Co-presence Frequency",
        title = "Top $(nrow(top_comms)) Present $K-Communities, thre=$(@sprintf("%.0e", thre))",
        markershape = :markershape,
        markerstrokecolor = :auto,
        lw = 2
    )
    display(plt1)
end

function summarize_size(df::DataFrame, taxa::Vector{String}, K::Int, thre::Float64, cal_counts::Dict)
    communities = collect(combinations(Symbol.(taxa), K));
    results = DataFrame(community = String[], freq_low = Float64[], freq_mid = Float64[], freq_high = Float64[], presence = Int[], peak = String[]);
    for comm in communities
        r = summarize_community(df, comm, cal_counts, thre)
        r !== nothing && push!(results, r)
    end
    results = transform(results, :peak => ByRow(p -> p == "mid" ? :solid : p == "high" ? :dash : :dot) => :linestyle);
    results = transform(results, :peak => ByRow(p -> p == "mid" ? :circle : p == "high" ? :cross : :diamond) => :markershape);
    return results
end

function summarize_system(df, taxa, K_ranges, thre, cal_counts, precut)
    all_counts = DataFrame();
    @showprogress for K in K_ranges
        results = summarize_size(df, taxa, K, thre, cal_counts);
        results.K = fill(K, nrow(results))
        results.presence_bin = ifelse.(results.presence .<=precut, "<=1", ">1")
        results.presence_bin = categorical(results.presence_bin; levels=["<=1", ">1"], ordered=true)

        counts_K = combine(groupby(results, [:K, :presence_bin, :peak]), nrow => :count)
        counts_K = combine(groupby(counts_K, [:K, :presence_bin])) do df
            total = sum(df.count)
            df.ratio = df.count ./ total
            df
        end

        counts_K.K = fill(string(K), nrow(counts_K)) # Store K as string for x-axis grouping
        append!(all_counts, counts_K)
    end
    all_counts.K = categorical(all_counts.K, ordered=true);
    all_counts.presence_bin = categorical(all_counts.presence_bin, levels=["<=1", ">1"], ordered=true);
    all_counts.peak = categorical(all_counts.peak, levels=["high", "mid", "low"], ordered=true);
    return all_counts
end

function show_all_counts(all_counts, thre)
    unique_K = sort!(unique(all_counts.K));
    presence_bins = ["<=1", ">1"];
    xvals = Dict{Tuple{String, String}, Float64}();
    offset = 0.0
    spacing_within = 0.4
    spacing_between = 1.0
    x_ticks = Float64[]
    x_tick_labels = String[]

    for k in unique_K
        for (i, bin) in enumerate(presence_bins)
            key = (string(k), bin)
            x = offset + (i-1) * spacing_within
            xvals[key] = x
            push!(x_ticks, x)
            push!(x_tick_labels, "$(k)-$bin")
        end
        offset += spacing_between
    end

    all_counts.x = [xvals[(string(row.K), row.presence_bin)] for row in eachrow(all_counts)];
    xtick_positions = [mean([xvals[(string(k), b)] for b in presence_bins]) for k in unique_K];
    xtick_labels = String.(unique_K);

    plt2 = plot(size=(400,250), grid=false, framestyle=:box);
    @df all_counts groupedbar!(
        plt2,
        :x,
        :ratio,
        group = :peak,
        bar_position = :stack,
        xlabel = "Size",
        ylabel = "Relative Count",
        guidefont=font(10),
        tickfont=font(8),
        legendfont=font(8),
        titlefont=font(12),
        bar_width = 0.3,
        legend= :outerright,
        xticks = (xtick_positions, xtick_labels),
        title = "Distributions of Peaks, threshold=$(@sprintf("%.0e", thre))"
    );

    display(plt2)
end


# -------------- Processing the dataset --------------
# --- Read metadata (calories, timestamps) ---
meta_keys = CSV.read("data/lifestyle/metadata/meta.keys", DataFrame; header=false)[!,2];
metadata = CSV.read("data/lifestyle/metadata/subjectA.gut.M.proc", DataFrame; header=false);
rename!(metadata, Symbol.(meta_keys));
calories = metadata[!, :nutrition_calorie];
cals_cut = quantile(calories[.!isnan.(calories)], [1/3, 2/3]);
timestamps = CSV.read("data/lifestyle/metadata/meta.time", DataFrame; header=false)[!,1];
time_calo = timestamps[findall(calories .> 0)];

# histogram(calories, bins=50, title="Calories per day", xlabel="Calories", ylabel="Frequency",
#     legend=false, grid=false, color=:blue, alpha=0.5)

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
microb_calo = subset!(microb_counts, :time => ByRow(in(time_calo)));
microb_calo.calories = [calories[findfirst(==(t), timestamps)] for t in microb_calo.time];
microb_calo = select(microb_calo, :time, :calories, Not([:time, :calories])); # reorder colums

# --- Aggregate OTU counts into taxonomy level ---
taxo_level = "Order" # choose from "Order", "Family", "Genus", "Species"
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
    for j in 3:ncol(rela_abun_calo)
        rela_abun_calo[i, j] /= total
    end
end;

# -------------- Anallyze the results --------------
# --- 1. top 10 presence communities ---
K = 4
final_rela_calo, taxa_candidates = filter_taxa(rela_abun_calo, 1e-3);
final_rela_calo, cal_counts = categorize_calories(final_rela_calo, cals_cut);
results_K4 = summarize_size(final_rela_calo, taxa_candidates, K, 1e-3, cal_counts);

histogram(results_K4.presence)
show_top_communities(results_K4, 20, K, 1e-3)

# ---- 2. Summarize peak distributions for K∈[2,3,4,5] ----
K_ranges = [2,3,4,5];

final_rela_calo, taxa_candidates = filter_taxa(rela_abun_calo, 1e-4);
final_rela_calo, cal_counts = categorize_calories(final_rela_calo, cals_cut);
all_counts = summarize_system(final_rela_calo, taxa_candidates, K_ranges, 1e-4, cal_counts, 10);

show_all_counts(all_counts, 1e-4)