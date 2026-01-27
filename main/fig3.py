import json
import numpy as np
import matplotlib.pyplot as plt

from helper import estimate_function_mean

target_lw = 0.5
params = {
    'font.size': 7,
    'axes.labelsize': 8,
    'xtick.labelsize': 7,
    'ytick.labelsize': 7,
    'axes.linewidth': target_lw,
    'xtick.major.width': target_lw,
    'ytick.major.width': target_lw,
}
plt.rcParams.update(params)


def make_grid_figure(json_path):
    with open(json_path, 'r') as f:
        data = json.load(f)
    
    S_list = [2, 4, 6, 8]
    
    def group_by_S(raw_list):
        grouped = {s: [] for s in S_list}
        for row in raw_list:
            if row['S'] in grouped:
                grouped[row['S']].append(row)
        return grouped

    dict_matr = group_by_S(data.get('data_matr', [])) 
    dict_init = group_by_S(data.get('data_init', []))
    
    fig, axes = plt.subplots(2, 4, figsize=(6.0, 3.0), dpi=300, 
                             sharex=True, sharey='row', 
                             constrained_layout=True)
    
    row_configs = [
        {
            'data': dict_init, 
            'color': 'darkorange', 
            'ylabel': r"Prob. Init. ($P_I$)",
            'ylim': (0.0, 1.0),
            'yticks': [0, 0.5, 1.0] 
        },
        {
            'data': dict_matr, 
            'color': 'dodgerblue', 
            'ylabel': r"Prob. Matr. ($P_M$)", 
            'ylim': (0.0, 0.55),
            'yticks': [0, 0.25, 0.5]
        }
    ]

    q_range = (1.0, 100.0)

    for row_idx, config in enumerate(row_configs):
        current_dict = config['data']
        color = config['color']
        
        for col_idx, S in enumerate(S_list):
            ax = axes[row_idx, col_idx]
            rows = current_dict.get(S, [])
            
            if rows:
                xs, p_mean, p_lower, p_upper = estimate_function_mean(rows, q_range)
                # mean_curve = np.exp(log_mean)
                # lower = np.exp(log_mean - 2 * log_std)
                # upper = np.exp(log_mean + 2 * log_std)
                
                ax.fill_between(xs, p_lower, p_upper, color=color, alpha=0.3, linewidth=0)
                ax.plot(xs, p_mean, color=color, linewidth=1.0)
            
            ax.set_xscale("log")
            ax.set_xlim(q_range)
            
            ax.tick_params(direction='out', top=False, right=False)
            ax.spines['top'].set_visible(False)
            ax.spines['right'].set_visible(False)
            ax.spines['left'].set_position(('outward', 5))
            ax.spines['bottom'].set_position(('outward', 5))
            
            
            if row_idx == 0:
                ax.set_title(f"$S={S}$", fontsize=8, pad=5, fontweight='bold')
            
            if col_idx == 0:
                ax.set_ylabel(config['ylabel'])
                if config['yticks']:
                    ax.set_yticks(config['yticks'])
                if config['ylim']:
                    ax.set_ylim(config['ylim'])


    fig.supxlabel(r"Total Energy Supply ($Q$)", fontsize=8)

    fig.align_ylabels(axes[:, 0])
    fig.savefig("figure3_grid_layout_v2.pdf")
    plt.show()

if __name__ == "__main__":
    make_grid_figure("data/output/fig3_data.json")