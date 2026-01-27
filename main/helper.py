import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.special import logit, expit

def estimate_function_mean(rows, x_range, n_points=300):
    xs = np.logspace(np.log10(x_range[0]), np.log10(x_range[1]), n_points)
    y_interp_list = []
    
    for row in rows:
        x = np.array(row['Qs'])
        y = np.array(row['y'])
        if len(x) < 2: continue 
        
        f = interp1d(x, y, kind='linear', bounds_error=False, fill_value=(y[0], y[-1]))
        y_interp_list.append(f(xs))
    
    if not y_interp_list: return xs, None, None

    y_stack = np.vstack(y_interp_list)
    
    # Clamping
    epsilon = 1e-3
    y_clamped = np.clip(y_stack, epsilon, 1.0 - epsilon)
    logit_y = logit(y_clamped) 
    
    mu_logit = np.mean(logit_y, axis=0)
    sigma_logit = np.std(logit_y, axis=0)
    
    
    p_median_curve = expit(mu_logit) 
    p_lower = expit(mu_logit - sigma_logit)
    p_upper = expit(mu_logit + sigma_logit)
    
    return xs, p_median_curve, p_lower, p_upper

