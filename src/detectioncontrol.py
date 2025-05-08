# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats

def condition_on_missingness(proteome):
    missing = proteome.isna().sum(axis = 0)
    missing = missing/proteome.shape[0]*100
    return missing.sort_values()

def calculate_optimum(data_in, alpha = .95):
    missingness = condition_on_missingness(data_in.proteome)

    # Parametric stats
    mu = np.mean(missingness)
    sigma = np.std(missingness, ddof=1)
    z = stats.norm.ppf(alpha)  # One-sided Z value
    upper_bound = mu + z * sigma

    upper_bound = np.percentile(missingness, alpha*100)

    # Plot empirical KDE with normal overlay
    plt.figure(figsize=(6, 4))
    sns.kdeplot(missingness, fill=True, label='Empirical KDE')

    # Plot normal PDF with same mean and std
    x_vals = np.linspace(min(missingness), max(missingness), 1000)
    normal_pdf = stats.norm.pdf(x_vals, mu, sigma)
    plt.plot(x_vals, normal_pdf, linestyle='--', label='Normal PDF', color='black')

    # Add vertical line at upper bound
    plt.axvline(upper_bound, color='red', linestyle='--', label=f'{int(alpha*100)}% Upper Bound = {upper_bound:.2f}%')

    # Labels and legend
    plt.xlabel('Sample Missingness (%)')
    plt.title('Missingness Distribution with Parametric CI')
    plt.legend()
    plt.tight_layout()
    plt.show()

    return upper_bound

def plot_retention(thresholds, response, optimal_threshold, file_out):
    plt.figure(figsize = (6,5))
    plt.plot(thresholds, response)
    plt.title(file_out)
    plt.xlabel('Sample Missingness Threshold %')
    plt.axvline(optimal_threshold, linestyle = '--', color = 'red')
    plt.savefig(file_out)
    plt.show()

def detection_control(data_in, optimal_threshold):
    missingness = condition_on_missingness(data_in.proteome)
    thresholds = np.linspace(0, 100, 21)
    recovered_counts = (missingness.values[None, :] <= thresholds[:, None]).sum(axis=1)
    plot_retention(thresholds, recovered_counts, optimal_threshold, os.path.join(data_in.outer_path, 'recovery_plot.pdf'))

    return data_in.proteome.loc[:,missingness<optimal_threshold]