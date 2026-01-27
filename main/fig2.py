import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.path as mpath

target_lw = 0.5
params = {
    'font.size': 7,
    'axes.labelsize': 8,
    'xtick.labelsize': 6,
    'ytick.labelsize': 6,
    'xtick.major.pad': 1,
    'ytick.major.pad': 1,
    'axes.linewidth': target_lw,
    'xtick.major.width': target_lw,
    'ytick.major.width': target_lw,
}
plt.rcParams.update(params)

class GeometryCalculator:
    def __init__(self, sigma, d, N0, Q_energy):
        self.sigma = np.array(sigma)
        self.d = np.array(d)
        self.N0 = np.array(N0)
        self.Q_energy = Q_energy
        
    def get_initialization_domain_polygon(self):
        vertices = [self.d]
        numerator = self.Q_energy - np.dot(self.N0, self.d)
        
        # 计算射线与 N0.s = Q 的交点
        intersections = []
        for i in range(2):
            v = self.sigma[:, i]
            denominator = np.dot(self.N0, v)
            if abs(denominator) > 1e-9:
                t = numerator / denominator
                if t > 0:
                    intersections.append(self.d + t * v)
        
        if len(intersections) == 2:
            vertices.extend(intersections)
            
        return np.array(vertices)

    def get_simplex_boundary(self):
        intercept_x = self.Q_energy / self.N0[0]
        intercept_y = self.Q_energy / self.N0[1]
        return np.array([[0, 0], [intercept_x, 0], [0, intercept_y]])

def calculate_maturation_ellipse(sigma, d, Q_threshold):
    sigma_mat = np.array(sigma)
    d_vec = np.array(d)
    inv_sigma = np.linalg.inv(sigma_mat)
    P = (inv_sigma + inv_sigma.T) / 2.0
    c = inv_sigma @ d_vec
    inv_P = np.linalg.inv(P)
    mu = 0.5 * inv_P @ c
    term_quadratic = mu.T @ P @ mu
    K = Q_threshold + term_quadratic
    
    if K <= 0: return None, None
    Sigma_plot = K * inv_P
    return mu, Sigma_plot

def get_ellipse_geometry(cov_matrix):
    vals, vecs = np.linalg.eigh(cov_matrix)
    order = vals.argsort()[::-1]
    vals = vals[order]
    vecs = vecs[:, order]
    
    theta = np.degrees(np.arctan2(vecs[1, 0], vecs[0, 0]))
    width, height = 2 * np.sqrt(vals)

    return width, height, theta


def plot_exact_geometry(ax, Q_val=10.0):
    sigma = [[1.45, 1.77], [0.45, 1.98]]
    d = [1.2, 1.8]
    N0 = [1.0, 1.2]
    
    calc = GeometryCalculator(sigma, d, N0, Q_val)
    
    gray_val = np.interp(Q_val, [7, 15], [0.97, 0.91])
    simplex_color = str(gray_val)
    simplex_verts = calc.get_simplex_boundary()
    simplex_path = mpath.Path(simplex_verts)
    simplex_patch = patches.PathPatch(simplex_path, facecolor=simplex_color, edgecolor='none', zorder=0)
    ax.add_patch(simplex_patch)


    cone_verts = calc.get_initialization_domain_polygon()
    poly_init = patches.Polygon(cone_verts, fc='none', ec='darkorange', lw=0, hatch='...', alpha=0.5, zorder=1)
    ax.add_patch(poly_init)
    
    if len(cone_verts) >= 3:
        ax.plot([cone_verts[0][0], cone_verts[1][0]], [cone_verts[0][1], cone_verts[1][1]], 
                color='#555555', lw=1, ls='-', zorder=2)
        ax.plot([cone_verts[0][0], cone_verts[2][0]], [cone_verts[0][1], cone_verts[2][1]], 
                color='#555555', lw=1, ls='-', zorder=2)

    intercepts = simplex_verts[1:]
    ax.plot([intercepts[0][0], intercepts[1][0]], 
            [intercepts[0][1], intercepts[1][1]], 
            color='darkorange', lw=1, zorder=3)

    matr_center, matr_cov = calculate_maturation_ellipse(sigma, d, Q_threshold=Q_val)
    
    if matr_center is not None:
        width, height, angle = get_ellipse_geometry(matr_cov)
        clip_cone = patches.Polygon(cone_verts, visible=False)
        ax.add_patch(clip_cone)
        clip_quadrant = patches.Rectangle((0,0), 100, 100, visible=False)
        ax.add_patch(clip_quadrant)

        ellipse_fill = patches.Ellipse(xy=matr_center, width=width, height=height, angle=angle,
                                       fc='dodgerblue', alpha=0.6, ec='none', zorder=3)
        ellipse_fill.set_clip_path(clip_cone)
        ax.add_patch(ellipse_fill)

        ellipse_edge = patches.Ellipse(xy=matr_center, width=width, height=height, angle=angle,
                                       fc='none', ec='dodgerblue', lw=1, zorder=4)
        ellipse_edge.set_clip_path(clip_quadrant)
        ax.add_patch(ellipse_edge)

    ax.scatter(d[0], d[1], s=10, c='gray', zorder=1.5)
    max_val = max(simplex_verts[1][0], simplex_verts[2][1])
    # ax.set_xlim(0, max_val * 1.05)
    # ax.set_ylim(0, max_val * 1.05)
    ax.set_xlim(0, 15.0)
    ax.set_ylim(0, 15.0)
    ax.set_aspect('equal')

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_position(('outward', 3))
    ax.spines['bottom'].set_position(('outward', 3))
    
    ax.set_xlabel(r"$s_1$")
    ax.set_ylabel(r"$s_2$")
    ax.set_title(f"$Q=${Q_val}", fontsize=7)

    ax.set_xticks([0, 15])
    ax.set_yticks([0, 15])


from matplotlib.lines import Line2D

def add_legend_A(ax):
    legend_elements = [
        (Line2D([0], [0], color='black', lw=0.5), "Nonnegative energy capture", 'black'),
        (Line2D([0], [0], color='#555555', lw=1), "Feasible steady-state biomass", '#555555'),
        (Line2D([0], [0], color='darkorange', lw=1), "Minimal energy capture", 'darkorange'),
        (Line2D([0], [0], color='dodgerblue', lw=1), "Steady-state energy capture", 'dodgerblue')
    ]
    
    handles = [h for h, l, c in legend_elements]
    labels = [l for h, l, c in legend_elements]
    colors = [c for h, l, c in legend_elements]

    leg = ax.legend(handles=handles, labels=labels, loc='upper right', frameon=False, 
                    fontsize=6, 
                    handlelength=1.0,
                    handletextpad=0.5,
                    borderaxespad=0,    
                    borderpad=2,        
                    labelspacing=0.4)
    
    for text, col in zip(leg.get_texts(), colors):
        text.set_color(col)

def add_legend_B(ax):
    handle_capture = patches.Patch(facecolor='#F0F0F0', label='Capture domain')
    handle_init = patches.Patch(facecolor='none', edgecolor='darkorange', 
                                hatch='...', alpha=0.5, label='Initialization domain')
    handle_matr = patches.Patch(facecolor='dodgerblue', alpha=0.6, label='Maturation domain')
    
    leg = ax.legend(handles=[handle_capture, handle_init, handle_matr], 
              loc='upper right', 
              fontsize=6, 
              frameon=False,
              handlelength=1.0, 
              handleheight=1.0,
              handletextpad=0.5,
              labelspacing=0.4,
              borderaxespad=0,
              borderpad=2)
    
    colors = [str(0.6), 'darkorange', 'dodgerblue']
    for text, col in zip(leg.get_texts(), colors):
        text.set_color(col)


def plot_result_curve(ax):
    try:
        data = np.load("data/output/fig2-example.npz")
        Q_range, prob_init, prob_matr = data["Q_range"], data["prob_init"], data["prob_matr"]
    except:
        Q_range = np.logspace(0, np.log10(250), 100)
        prob_matr = 0.25 * (1 - np.exp(-Q_range/20))
        prob_init = prob_matr * 1.3 + 0.02

    ax.plot(Q_range, prob_init, color='darkorange', label="Prob. Initialization", 
            lw=1.2, ls=':', zorder=2)
    ax.plot(Q_range, prob_matr, color='dodgerblue', label="Prob. Maturation", 
            lw=1.2, ls='-', zorder=3)

    target_Qs = [7, 15]
    labels = ["A", "B"]
    y_label_level = 0.15
    
    for q, lab in zip(target_Qs, labels):
        y_init = np.interp(q, Q_range, prob_init)
        y_matr = np.interp(q, Q_range, prob_matr)

        ax.scatter([q, q], [y_init, y_matr], color=['darkorange', 'dodgerblue'], 
                   s=15, zorder=4, edgecolors='white', linewidths=0.5)
        
        # ax.text(q, y_label_level, f"({lab})", ha='center', va='center',
        #         fontsize=8, fontweight='bold', color='#333333',
        #         bbox=dict(facecolor='white', edgecolor='none', alpha=0.7, pad=1), # 增加微型背景防止干扰
        #         zorder=5)


    ax.legend(loc='upper left', fontsize=7, frameon=False)
    ax.set_xscale('log')
    ax.set_xlim(1, 250)
    ax.set_ylim(0, 0.4)
    ax.set_yticks([0.15, 0.3])
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.set_xlabel(r"Total Energy Supply $(Q)$", labelpad=1) 
    ax.set_ylabel(r"Probability", labelpad=1)


fig = plt.figure(figsize=(4.5, 4), dpi=300, layout='constrained')
gs = fig.add_gridspec(2, 2, height_ratios=[1.2, 1])

ax1 = fig.add_subplot(gs[0, 0])
plot_exact_geometry(ax1, Q_val=7)
# ax1.set_title("A", loc='left', fontweight='bold') # Add panel label
add_legend_A(ax1)

ax2 = fig.add_subplot(gs[0, 1])
plot_exact_geometry(ax2, Q_val=15)
# ax2.set_title("B", loc='left', fontweight='bold') # Add panel label
add_legend_B(ax2)

ax3 = fig.add_subplot(gs[1, :])
plot_result_curve(ax3)

plt.savefig("figures/fig2-example.pdf", dpi=300)
plt.show()