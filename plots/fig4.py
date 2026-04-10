import json
import matplotlib.pyplot as plt

from helper import estimate_function_mean, estimate_curve_mean_std


target_lw = 0.5
params = {
    "font.size": 7,
    "axes.labelsize": 8,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "axes.linewidth": target_lw,
    "xtick.major.width": target_lw,
    "ytick.major.width": target_lw,
}
plt.rcParams.update(params)


def _style_axis(ax):
    ax.set_xscale("log")
    ax.tick_params(direction="out", top=False, right=False)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_position(("outward", 5))
    ax.spines["bottom"].set_position(("outward", 5))


def _pm_size_colors(groups):
    palette = [
        "#7D88C7",  # blue-violet
        "#6E97D6",  # periwinkle blue
        "#69C2D0",  # blue-cyan
        "#67CDBA",  # cyan-teal
        "#6FAFDB",  # airy blue
        "dodgerblue",
    ]
    sizes = [group["comm_size"] for group in groups]
    colors = {}
    for idx, size in enumerate(sizes):
        colors[size] = palette[min(idx, len(palette) - 1)]
    if sizes:
        colors[max(sizes)] = "dodgerblue"
    return colors


def _load_fig4_data(json_path):
    with open(json_path, "r") as f:
        return json.load(f)


def _prepare_fig4_data(data):
    meta = data["meta"]
    pm_groups = sorted(data["pm_groups"], key=lambda row: row["comm_size"])
    pm_rows = data["pm_rows"]
    richness_rows = data["richness_rows"]
    pm_colors = _pm_size_colors(pm_groups)

    rows_by_size = {group["comm_size"]: [] for group in pm_groups}
    for row in pm_rows:
        rows_by_size[row["comm_size"]].append(row)

    return {
        "S": meta["S"],
        "x_range": tuple(meta["x_range"]),
        "pm_ymax": meta["pm_ymax"],
        "richness_ymax": meta["richness_ymax"],
        "pm_groups": pm_groups,
        "rows_by_size": rows_by_size,
        "richness_rows": richness_rows,
        "pm_colors": pm_colors,
    }


def _resolve_by_S(default, S, overrides):
    if overrides is None:
        return default
    if isinstance(overrides, dict):
        return overrides.get(S, default)
    return overrides


def _plot_pm_panel(ax, plot_data, pm_ymax=None, pm_yticks=None, show_legend=True):
    x_range = plot_data["x_range"]
    pm_groups = plot_data["pm_groups"]
    rows_by_size = plot_data["rows_by_size"]
    pm_colors = plot_data["pm_colors"]

    for group in pm_groups:
        rows = rows_by_size[group["comm_size"]]
        xs, mean, lower, upper = estimate_function_mean(rows, x_range)
        if mean is None:
            continue

        color = pm_colors[group["comm_size"]]
        ax.fill_between(xs, lower, upper, color=color, alpha=0.20, linewidth=0)
        ax.plot(
            xs,
            mean,
            color=color,
            linewidth=1.2,
            label=f"{group['comm_size']}",
        )

    if pm_ymax is None:
        pm_ymax = plot_data["pm_ymax"]
    if pm_yticks is None:
        pm_yticks = [0.0, 0.05, 0.10]

    _style_axis(ax)
    ax.set_xlim(x_range)
    ax.set_ylim((0.0, pm_ymax))
    ax.set_yticks(pm_yticks)
    ax.set_ylabel(r"Prob. Matr. ($P_M$)")

    handles, labels = ax.get_legend_handles_labels()
    if show_legend and handles:
        legend = ax.legend(
            title="Comm. Size",
            title_fontsize=6,
            frameon=False,
            fontsize=5,
            loc="upper right",
            bbox_to_anchor=(0.93, 1.00),
            handlelength=1.75,
            handletextpad=0.75,
            labelspacing=0.1,
            borderpad=0.25,
        )
        legend._legend_box.align = "left"


def _plot_richness_panel(ax, plot_data, richness_ymax=None, richness_yticks=None):
    x_range = plot_data["x_range"]
    richness_rows = plot_data["richness_rows"]

    for row in richness_rows:
        ax.plot(
            row["Qs"],
            row["y"],
            color="#88A5BD",
            alpha=0.18,
            linewidth=0.5,
        )

    xs, mean, lower, upper = estimate_curve_mean_std(richness_rows, x_range)
    if mean is not None:
        ax.fill_between(xs, lower, upper, color="#3F78A8", alpha=0.10, linewidth=0)
        ax.plot(xs, mean, color="#E4EFF8", linewidth=3.8)
        ax.plot(xs, mean, color="#2E5F8A", linewidth=1.8)

    if richness_ymax is None:
        richness_ymax = plot_data["richness_ymax"]
    if richness_yticks is None:
        richness_yticks = [0, 2, 4]

    _style_axis(ax)
    ax.set_xlim(x_range)
    ax.set_ylim((0.0, richness_ymax))
    ax.set_yticks(richness_yticks)
    ax.set_ylabel(r"Ave. Diversity ($\langle D \rangle$)")


def _plot_fig4_row(
    axes,
    plot_data,
    pm_ymax=None,
    richness_ymax=None,
    pm_yticks=None,
    richness_yticks=None,
    show_legend=True,
):
    _plot_pm_panel(
        axes[0],
        plot_data,
        pm_ymax=pm_ymax,
        pm_yticks=pm_yticks,
        show_legend=show_legend,
    )
    _plot_richness_panel(
        axes[1],
        plot_data,
        richness_ymax=richness_ymax,
        richness_yticks=richness_yticks,
    )


def make_fig4(
    json_path,
    outpath="figures/Figure4.pdf",
    pm_ymax=None,
    richness_ymax=None,
    pm_yticks=None,
    richness_yticks=None,
):
    plot_data = _prepare_fig4_data(_load_fig4_data(json_path))

    fig, axes = plt.subplots(
        1,
        2,
        figsize=(6.0, 2.1),
        dpi=300,
        constrained_layout=True,
    )

    _plot_fig4_row(
        axes,
        plot_data,
        pm_ymax=pm_ymax,
        richness_ymax=richness_ymax,
        pm_yticks=pm_yticks,
        richness_yticks=richness_yticks,
        show_legend=True,
    )

    fig.supxlabel(r"Total Energy Supply ($Q$)", fontsize=8)
    fig.align_ylabels(axes)
    fig.savefig(outpath)
    plt.show()


def make_fig4_si(
    json_paths=None,
    S_list=(2, 4, 8),
    outpath="figures/Figure4_SI.pdf",
    pm_ymax_by_S={2: 0.60, 4: 0.22, 8: 0.05},
    pm_yticks_by_S={2: [0.0, 0.3, 0.6], 4: [0.0, 0.1, 0.2], 8: [0.0, 0.025, 0.05]},
    richness_ymax_by_S=None,
    richness_yticks_by_S=None,
):
    if json_paths is None:
        json_paths = {S: f"data/output/fig4_data-S{S}.json" for S in S_list}
    elif not isinstance(json_paths, dict):
        json_paths = dict(zip(S_list, json_paths))

    plot_data_by_S = {}
    for S in S_list:
        plot_data = _prepare_fig4_data(_load_fig4_data(json_paths[S]))
        plot_data_by_S[S] = plot_data

    height_ratios = []
    for _ in S_list:
        height_ratios.extend([0.16, 1.0])

    fig = plt.figure(
        figsize=(6.0, 2.15 * len(S_list)),
        dpi=300,
        constrained_layout=True,
    )
    gs = fig.add_gridspec(
        2 * len(S_list),
        2,
        height_ratios=height_ratios,
    )

    axes = []
    shared_pm_ax = None
    shared_richness_ax = None

    for row_idx, S in enumerate(S_list):
        title_ax = fig.add_subplot(gs[2 * row_idx, :])
        title_ax.axis("off")
        title_ax.text(
            0.5,
            0.0,
            f"$S={plot_data_by_S[S]['S']}$",
            ha="center",
            va="bottom",
            fontsize=8,
            fontweight="bold",
        )

        pm_ax = fig.add_subplot(
            gs[2 * row_idx + 1, 0],
            sharex=shared_pm_ax,
        )
        richness_ax = fig.add_subplot(
            gs[2 * row_idx + 1, 1],
            sharex=shared_richness_ax,
        )
        if shared_pm_ax is None:
            shared_pm_ax = pm_ax
        if shared_richness_ax is None:
            shared_richness_ax = richness_ax
        axes.append((pm_ax, richness_ax))

        plot_data = plot_data_by_S[S]
        _plot_fig4_row(
            axes[row_idx],
            plot_data,
            pm_ymax=_resolve_by_S(plot_data["pm_ymax"], S, pm_ymax_by_S),
            richness_ymax=_resolve_by_S(plot_data["richness_ymax"], S, richness_ymax_by_S),
            pm_yticks=_resolve_by_S([0.0, 0.05, 0.10], S, pm_yticks_by_S),
            richness_yticks=_resolve_by_S([0, 2, 4], S, richness_yticks_by_S),
            show_legend=True,
        )

    fig.supxlabel(r"Total Energy Supply ($Q$)", fontsize=8)
    fig.align_ylabels([row[0] for row in axes])
    fig.align_ylabels([row[1] for row in axes])
    fig.savefig(outpath)
    plt.show()


if __name__ == "__main__":
    # make_fig4("data/output/fig4_data.json")
    make_fig4_si()
