import json
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from helper import estimate_function_mean


target_lw = 0.5
draw_insect = False
params = {
    'font.size': 7,
    'axes.labelsize': 8,        
    'xtick.labelsize': 6,      
    'ytick.labelsize': 6,       
    'axes.titlesize': 8,
    'axes.linewidth': target_lw,
    'xtick.major.width': target_lw,
    'ytick.major.width': target_lw,
}
plt.rcParams.update(params)


def plot_embedded_figure(json_path, q_range=(1.0, 1000.0)):
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    curve_groups = {}
    for row in data:
        k = row['k']
        if k not in curve_groups: curve_groups[k] = []
        curve_groups[k].append(row)

    target_ks = [0.0, 0.2, 0.4]
    ncols = len(target_ks)
    
    fig, axes = plt.subplots(1, ncols, figsize=(6.0, 2.0), dpi=300, constrained_layout=True, sharey=True)
    
    curve_colors = [cm.Blues(x) for x in np.linspace(0.8, 0.4, ncols)]
    cf_handle = None

    for col, k_val in enumerate(target_ks):
        ax = axes[col]
        
        rows = curve_groups.get(k_val, [])
        cl = curve_colors[col]
        if rows:
            xs, p_mean, p_lower, p_upper = estimate_function_mean(rows, q_range)
            ax.fill_between(xs, p_lower, p_upper, color=cl, alpha=0.2, linewidth=0)
            ax.plot(xs, p_mean, color=cl, linewidth=1.5)

        ax.set_xscale("log")
        ax.set_xlim(q_range)
        ax.set_title(f"$k={k_val}$", fontweight='bold', pad=10)
        
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_position(('outward', 5))
        ax.spines['bottom'].set_position(('outward', 5))
        ax.tick_params(direction='out', top=False, right=False)

        if col == 0:
            ax.set_ylabel(r"Prob. Maturation ($P_M$)")
            ax.set_ylim(0.0, 0.75)
            ax.set_yticks([0.0, 0.2, 0.4, 0.6])
        else:
            ax.tick_params(labelleft=False)


        if draw_insect:
            ax_ins = ax.inset_axes([0.6, 0.63, 0.35, 0.35])
            
            filename = f"data/output/Q_contour_k={k_val}.npz"
            try:
                d = np.load(filename)
                s1, s2, Q = d["s1"], d["s2"], d["Q"]
                S1, S2 = np.meshgrid(s1, s2, indexing="ij")
            except:
                s1 = np.linspace(0, 6, 100); s2 = np.linspace(0, 6, 100)
                S1, S2 = np.meshgrid(s1, s2, indexing="ij")
                Q = np.sqrt((S1-1)**2 + (S2-1)**2) * (1 if k_val==0.0 else 0.8)

            cf = ax_ins.contourf(S1, S2, Q, levels=200, cmap="GnBu", vmin=0.0, vmax=10.0)
            if col == 2: cf_handle = cf

            ax_ins.set_aspect("equal")
            ax_ins.set_xlim(0.8, 8.2)
            ax_ins.set_ylim(0.8, 8.2)
            
            ax_ins.set_xticks([])
            ax_ins.set_yticks([]) 
            ax_ins.tick_params(bottom=False, left=False)
            ax_ins.spines['top'].set_visible(False)
            ax_ins.spines['right'].set_visible(False)
            ax_ins.spines['bottom'].set_visible(True)
            ax_ins.spines['left'].set_visible(True)
        
            ax_ins.patch.set_alpha(0.8)
        
            ax_ins.set_xlabel(r"$\boldsymbol{s}_1$", fontsize=6, labelpad=1)
            ax_ins.set_ylabel(r"$\boldsymbol{s}_2$", fontsize=6, labelpad=1)
            
            for spine in ax_ins.spines.values():
                spine.set_linewidth(0.3)
            
            # cbar = fig.colorbar(cf_handle, ax=axes, location='right', 
                                # shrink=0.4, aspect=15, pad=0.02)

    cax = axes[2].inset_axes([1.08, 0.63, 0.03, 0.3])
    cbar = fig.colorbar(cf_handle, cax=cax)

    cbar.set_ticks([0, 5, 10])
    cbar.set_label(r"$Q(\boldsymbol{s})$", rotation=270, labelpad=5, fontsize=6)
    cbar.outline.set_linewidth(target_lw)
    cbar.ax.tick_params(width=target_lw, labelsize=5)

    fig.set_constrained_layout_pads(w_pad=0.07, h_pad=0.07, hspace=0, wspace=0.1)
    fig.supxlabel(r"Total Energy Supply ($Q$)", fontsize=8)
    
    plt.savefig("figures/saturation-S4.pdf")
    plt.show()

# S = 4
plot_embedded_figure("data/output/k_curves_data-S4.json")

# S = 2
# plot_embedded_figure("data/output/k_curves_data-S2.json")
