# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
import pandas as pd
import numpy as np

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from carrier import Carrier
import pickle
import argparse

def condition_on_missingness(carrier):
    missing = carrier.proteome.isna().sum(axis=0)
    missing = missing / carrier.proteome.shape[0] * 100
    carrier.missingness = missing.sort_values()
    return carrier.missingness

def calculate_optimum_threshold(carrier, alpha=0.99):
    missingness = condition_on_missingness(carrier)

    # Parametric stats
    mu = np.mean(missingness)
    sigma = np.std(missingness, ddof=1)
    z = stats.norm.ppf(alpha)
    upper_bound = mu + z * sigma

    # Plot KDE with normal PDF overlay
    PDF_plot = plt.figure(figsize=(6, 5))
    x_plot = missingness.dropna()  
    try:
        sns.kdeplot(x_plot, fill=True, label='Empirical KDE')   # seaborn ≥ 0.11
    except TypeError:
        sns.kdeplot(x_plot, shade=True, label='Empirical KDE')  # seaborn < 0.11
    x_vals = np.linspace(float(x_plot.min()), float(x_plot.max()), 1000) if len(x_plot) else np.linspace(0, 100, 1000)
    normal_pdf = stats.norm.pdf(x_vals, mu, sigma)
    plt.plot(x_vals, normal_pdf, linestyle='--', label='Normal PDF', color='black')
    plt.axvline(upper_bound, color='red', linestyle='--', label=f'{int(alpha * 100)}% Upper Bound = {upper_bound:.2f}%')
    plt.xlabel('Sample Missingness (%)')
    plt.title(f'{carrier.projectname}\nMissingness PDF with Parametric CI (alpha = {alpha})')
    plt.legend()
    plt.tight_layout()
    # plt.savefig(os.path.join(carrier.outerpath, carrier.projectname + '_' + str(alpha) + '_recovery_PDF.pdf'))
    carrier.PDF_plot = PDF_plot

    # Plot pseudo CDF
    thresholds = np.linspace(0, 100, 21)
    recovered_counts = (missingness.values[None, :] <= thresholds[:, None]).sum(axis=1)
    CDF_plot = plt.figure(figsize=(6, 5))
    plt.plot(thresholds, recovered_counts)
    plt.title(f'{carrier.projectname}\nPseudo CDF using sample retention (alpha = {alpha})')
    plt.xlabel('Sample Missingness Threshold %')
    plt.axvline(upper_bound, linestyle='--', color='red')
    plt.tight_layout()
    #plt.savefig(os.path.join(carrier.outerpath, carrier.projectname + '_' + str(alpha) + '_recovery_pseudoCDF.pdf'))
    carrier.CDF_plot = CDF_plot
    
    return upper_bound

def detection_control(carrier, optimal_threshold):
    carrier.proteome = carrier.proteome.loc[:, carrier.missingness < optimal_threshold]
    carrier.status = f'dcontrol_{int(optimal_threshold)}'
    print(f"Dropped {(carrier.missingness >= optimal_threshold).sum()} samples (kept {(carrier.missingness < optimal_threshold).sum()})")
    return carrier

def main(alpha):
    ALPHA = alpha

    with open(os.path.join(os.path.dirname(__file__), "../output", "shogoki.pkl"), "rb") as f:
        nigoki = pickle.load(f)

    print("3. Running detection control...")
    optimum = calculate_optimum_threshold(nigoki, alpha=ALPHA)
    detection_control(nigoki, optimum)

    outdir = nigoki.outerpath
    if not isinstance(outdir, str) or not os.path.isdir(outdir):
        outdir = os.path.abspath(os.path.join(os.path.dirname(__file__), "../output"))
    os.makedirs(outdir, exist_ok=True)

    # write to disk
    with open(os.path.join(outdir, 'nigoki.pkl'), 'wb') as f:
        pickle.dump(nigoki, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "-a", type=float, default=0.99, required=False,
        help="Alpha/significance parameter for detection control (default: 0.99)."
    )
    args = parser.parse_args()
    main(alpha=args.a)
