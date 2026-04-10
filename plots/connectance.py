import json
import matplotlib.pyplot as plt

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from plotting.helper import estimate_function_mean

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

C_VALUES = [1.0, 0.8, 0.6, 0.4, 0.2]
COLORS = ["#1f4e79", "#2e75b6", "#5ba3d9", "#9dc3e6", "#bdd7ee"]

json_path = "data/connectance/collected.json"
with open(json_path) as f:
    all_rows = json.load(f)

# group by connectance
by_c = {c: [] for c in C_VALUES}
for row in all_rows:
    c = round(row["c"], 2)
    if c in by_c:
        by_c[c].append({"Qs": row["Qs"], "y": row["pm"]})

q_range = (7.5, 100.0)

fig, ax = plt.subplots(1, 1, figsize=(3.0, 2.5), dpi=300, constrained_layout=True)

for (c, color) in zip(C_VALUES, COLORS):
    xs, mean, lower, upper = estimate_function_mean(by_c[c], q_range)
    if mean is None:
        continue
    ax.fill_between(xs, lower, upper, color=color, alpha=0.25, linewidth=0)
    ax.plot(xs, mean, color=color, linewidth=1.2, label=f"c={c}")

ax.set_xscale("log")
ax.set_xlim(q_range)
ax.tick_params(direction="out", top=False, right=False)
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["left"].set_position(("outward", 5))
ax.spines["bottom"].set_position(("outward", 5))
ax.set_ylabel(r"Prob. Matr. ($P_M$)")
ax.set_xlabel(r"Total Energy Supply ($Q$)")
ax.legend(frameon=False, fontsize=6, loc="upper right")

fig.savefig("figures/connectance.pdf")
plt.show()
