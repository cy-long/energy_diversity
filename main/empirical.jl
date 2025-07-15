using CSV, DataFrames, StatsBase, Statistics
using Plots, Random, Combinatorics

# ---- Read metadata (calories, timestamps) ----
meta_keys = CSV.read("data/lifestyle/metadata/meta.keys", DataFrame; header=false)[!,2];
metadata = CSV.read("data/lifestyle/metadata/subjectA.gut.M.proc", DataFrame; header=false);
rename!(metadata, Symbol.(meta_keys));
calories = metadata[!, :nutrition_calorie];

timestamps = CSV.read("data/lifestyle/metadata/meta.time", DataFrame; header=false)[!,1];
time_calo = timestamps[findall(calories .> 0)];

histogram(calories, bins=50, title="Calories per day", xlabel="Calories", ylabel="Frequency",
    legend=false, grid=false, color=:blue, alpha=0.5);


# ---- Read taxanomy ----
microb_keys = CSV.read("data/lifestyle/microbiota/subjectA.gut.keys", DataFrame; header=true);
taxo = split.(microb_keys.Taxo, ';');

microb_keys.Order   = replace.(getindex.(taxo, 4), r"^o__" => "");
microb_keys.Family  = replace.(getindex.(taxo, 5), r"^f__" => "");
microb_keys.Genus   = replace.(getindex.(taxo, 6), r"^g__" => "");
microb_keys.Species = replace.(getindex.(taxo, 7), r"^s__" => "");
select!(microb_keys, Not(:Taxo, :Source));


# ---- Read sequencing time series ----
microb_counts = CSV.read("data/lifestyle/microbiota/subjectA.gut.micro", DataFrame; header=false);

# rename columns (time and OTUs)
rename!(microb_counts, vcat(:time, Symbol.(microb_keys.OTU)));
microb_counts.time = Int.(microb_counts.time);

# filter rows (days) with calories records and combine them
microb_calo = subset!(microb_counts, :time => ByRow(in(time_calo)));
microb_calo.calories = [calories[findfirst(==(t), timestamps)] for t in microb_calo.time];
microb_calo = select(microb_calo, :time, :calories, Not([:time, :calories])); # reorder colums


# ---- Aggregate OTU counts into taxonomy level ----
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

# convert into relative abundance
rela_abun_calo = copy(taxa_calo);
for i in 1:nrow(rela_abun_calo)
    total = sum(rela_abun_calo[i, 3:end])
    for j in 3:ncol(rela_abun_calo)
        rela_abun_calo[i, j] /= total
    end
end;


# ---- Build presence-absence table; Select candidates ----
thre = 1e-3 # threshold for relative abundance

# exclude taxa that are always below or above the relative abundance threshold
remove_taxa = String[];
for col in names(rela_abun_calo)[3:end]
    vals = rela_abun_calo[!, col]
    if col == "" || all(<=(thre), vals) || all(>=(thre), vals)
        push!(remove_taxa, col)
    end
end

final_rela_calo = select(rela_abun_calo, Not(remove_taxa));
taxa_candidates = names(final_rela_calo)[3:end];
print("Final candidates under threshold $thre: ", taxa_candidates, "\n");

# ---- Analyze co-presence frequency vs. calories ----
K = 3 # combination size
cals_q = [2000.0, 2500.0]; # cals_q = quantile(calories[.!isnan.(calories)], [1/3, 2/3]);

is_low = final_rela_calo.calories .< cals_q[1]; n_low = sum(is_low)
is_mid = (final_rela_calo.calories .>= cals_q[1]) .& (final_rela_calo.calories .<= cals_q[2]); n_mid = sum(is_mid)
is_high = final_rela_calo.calories .> cals_q[2]; n_high = sum(is_high)
print("Days in low / mid / high calories: $n_low / $n_mid / $n_high\n");

# generate K-taxa communities
communities = collect(combinations(Symbol.(taxa_candidates), K));
results = DataFrame(community = String[], calorie_level = String[], frequency = Float64[]);

# calculate co-presence frequency for each community
for comm in communities
    comm_calo = select(final_rela_calo, [:calories, comm...])

    # count rows where all taxa > threshold (co-presence)
    copresent_low = sum(all(row .> thre) for row in eachrow(Matrix(comm_calo[is_low, 2:end])))
    copresent_mid = sum(all(row .> thre) for row in eachrow(Matrix(comm_calo[is_mid, 2:end])))
    copresent_high = sum(all(row .> thre) for row in eachrow(Matrix(comm_calo[is_high, 2:end])))

    # skip if the community never co-presents
    if copresent_low == 0 && copresent_mid == 0 && copresent_high == 0
        continue
    end

    # add to results
    comm_name = join(string.(comm), "+")
    push!(results, (comm_name, "low", copresent_low / n_low))
    push!(results, (comm_name, "mid", copresent_mid / n_mid))
    push!(results, (comm_name, "high", copresent_high / n_high))
end


# Plot: each combination is one line, x-axis is low/mid/high
plt1 = plot(legend=false); cls = ["low", "mid", "high"];
for comm in groupby(results, :community)
    plot!(cls, comm.frequency, label = false, marker = :circle, alpha = 0.5)
end

xlabel!("Calorie Level"); ylabel!("Co-presence Frequency");
display(plt1)

# Next to do: change K (size of combinations) and ><thre (initiation vs. realization)

# # --- Histogram: which calorie level lead to peak co-presence frequency ---
# max_levels = Int[];
# for grp in groupby(results, :combo)
#     push!(max_levels, argmax(grp.freq))
# end
# level_labels = ["low", "mid", "high"]
# counts = countmap(max_levels)

# # Bar plot on same axis
# bar(level_labels, [counts[i] for i in 1:3],
#     title = "Histogram of max coexistence level",
#     xlabel = "Calories range",
#     ylabel = "Number of combinations",
#     legend = false)

# ---- Let's do sth else: increasing the bins -----
# for combo in taxa_combos
#     sub = select(rela_abun_calo, [:calories, combo...])
#     cals = sub[:, 1]
#     coex_mask = [all(row .> thre) for row in eachrow(Matrix(sub[:, 2:end]))]

#     cals_coex = cals[coex_mask]
#     if isempty(cals_coex)
#         continue
#     end

#     histogram(cals_coex, bins=20, legend=false,
#               title=join(string.(combo), "+"),
#               xlabel="Calories", ylabel="Coexistence count")
#     display(current())  # or display(plot!)

#     break  # 🛑 stop after first plot
# end

# # ---- 1. maximal numbers of taxa at each different calories level ---
# thre = 1e-4
# taxa_data = Tables.matrix(rela_abun_calo[:, 3:end]);
# present = taxa_data .> thre;
# coexist_counts = sum(present, dims=2);
# coexist_counts = vec(coexist_counts);
# daily_cals = rela_abun_calo.calories;