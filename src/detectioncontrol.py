# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

def condition_on_missingness(data_in):

    missing = data_in.proteome.isna().sum(axis = 0)
    missing = missing/data_in.proteome.shape[0]*100
    data_in.missingness = missing.sort_values()
    return data_in.missingness

def calculate_optimum_threshold(data_in, alpha = .95):

    missingness = condition_on_missingness(data_in)

    # Parametric stats
    mu = np.mean(missingness)
    sigma = np.std(missingness, ddof=1)
    z = stats.norm.ppf(alpha)  # One-sided Z value
    upper_bound = mu + z * sigma

    # Plot empirical KDE with normal overlay
    plt.figure(figsize=(6,5))
    sns.kdeplot(missingness, fill=True, label='Empirical KDE')

    # Plot normal PDF with same mean and std
    x_vals = np.linspace(min(missingness), max(missingness), 1000)
    normal_pdf = stats.norm.pdf(x_vals, mu, sigma)
    plt.plot(x_vals, normal_pdf, linestyle='--', label='Normal PDF', color='black')
    # Add vertical line at upper bound
    plt.axvline(upper_bound, color='red', linestyle='--', label=f'{int(alpha*100)}% Upper Bound = {upper_bound:.2f}%')
    # Labels and legend
    plt.xlabel('Sample Missingness (%)')
    plt.title(f'{data_in.projectname}\nMissingness PDF with Parametric CI (alpha = {alpha})')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(data_in.outer_path, data_in.projectname + str(alpha) + '_recovery_PDF.pdf'))
    plt.show()
    
    # Plot pseudo CDF using cumulative retention of samples (not %)
    thresholds = np.linspace(0, 100, 21)
    recovered_counts = (missingness.values[None, :] <= thresholds[:, None]).sum(axis=1)
    plt.figure(figsize = (6,5))
    plt.plot(thresholds, recovered_counts)
    plt.title(f'{data_in.projectname}\nPseudo CDF using sample retention (alpha = {alpha})')
    plt.xlabel('Sample Missingness Threshold %')
    plt.axvline(upper_bound, linestyle = '--', color = 'red')
    plt.savefig(os.path.join(data_in.outer_path, data_in.projectname + str(alpha) + '_recovery_pseudoCDF.pdf'))
    plt.show()

    return upper_bound

def detection_control(data_in, optimal_threshold):

    data_in.proteome = data_in.proteome.loc[:,data_in.missingness<optimal_threshold]
    data_in.status = f'dcontrol_{int(optimal_threshold)}'

    return data_in