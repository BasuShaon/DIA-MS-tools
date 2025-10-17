# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from carrier import Carrier
import pickle
import argparse

def condition_on_missingness(carrier):
    """ 
    Computes the percentage of missing values per sample in the proteome and stores the sorted result.

    Parameters
    ----------
    carrier : object
        Data object that includes a 'proteome' DataFrame with precursor intensity data.

    Returns
    -------
    pd.Series
        Sorted percentage of missing values per sample, stored in `carrier.missingness`.
    """
    missing = carrier.proteome.isna().sum(axis=0)
    missing = missing / carrier.proteome.shape[0] * 100
    carrier.missingness = missing.sort_values()

    return carrier.missingness

def calculate_optimum_threshold(carrier, alpha=0.99):
    """ 
    Calculates the optimal missingness threshold using a parametric confidence interval and visualizes results.

    Parameters
    ----------
    carrier : object
        Data object that includes 'proteome', 'projectname', and 'outerpath'.
    alpha : float, optional
        Confidence level for calculating the upper bound threshold (default is 0.95).

    Returns
    -------
    upper_bound : float
        Optimal missingness threshold based on the upper bound of the normal distribution.
    """
    missingness = condition_on_missingness(carrier)

    # Parametric stats
    mu = np.mean(missingness)
    sigma = np.std(missingness, ddof=1)
    z = stats.norm.ppf(alpha)
    upper_bound = mu + z * sigma

    # Plot KDE with normal PDF overlay
    plt.figure(figsize=(6, 5))
    sns.kdeplot(missingness, fill=True, label='Empirical KDE')
    x_vals = np.linspace(min(missingness), max(missingness), 1000)
    normal_pdf = stats.norm.pdf(x_vals, mu, sigma)
    plt.plot(x_vals, normal_pdf, linestyle='--', label='Normal PDF', color='black')
    plt.axvline(upper_bound, color='red', linestyle='--', label=f'{int(alpha * 100)}% Upper Bound = {upper_bound:.2f}%')
    plt.xlabel('Sample Missingness (%)')
    plt.title(f'{carrier.projectname}\nMissingness PDF with Parametric CI (alpha = {alpha})')
    plt.legend()
    plt.tight_layout()
    plt.savefig(os.path.join(carrier.outerpath, carrier.projectname + '_' + str(alpha) + '_recovery_PDF.pdf'))
    # plt.show()

    # Plot pseudo CDF
    thresholds = np.linspace(0, 100, 21)
    recovered_counts = (missingness.values[None, :] <= thresholds[:, None]).sum(axis=1)
    plt.figure(figsize=(6, 5))
    plt.plot(thresholds, recovered_counts)
    plt.title(f'{carrier.projectname}\nPseudo CDF using sample retention (alpha = {alpha})')
    plt.xlabel('Sample Missingness Threshold %')
    plt.axvline(upper_bound, linestyle='--', color='red')
    plt.savefig(os.path.join(carrier.outerpath, carrier.projectname + '_' + str(alpha) + '_recovery_pseudoCDF.pdf'))
    # plt.show()

    return upper_bound


def detection_control(carrier, optimal_threshold):
    """ 
    Filters proteome data by removing samples above the missingness threshold.

    Parameters
    ----------
    carrier : object
        Data object containing the proteome and computed missingness.
    optimal_threshold : float
        Threshold for allowed percentage of missing values per sample.

    Returns
    -------
    object
        Modified carrier object with filtered `proteome` and updated `status`.
    """

    carrier.proteome = carrier.proteome.loc[:, carrier.missingness < optimal_threshold]
    carrier.status = f'dcontrol_{int(optimal_threshold)}'

    return carrier

def main(alpha):
    # Configuration
    ALPHA = alpha

    with open(os.path.join(os.path.dirname(__file__), "../output", "shogoki.pkl"), "rb") as f:
        nigoki = pickle.load(f)

    print("3. Running detection control...")
    optimum = calculate_optimum_threshold(nigoki, alpha=ALPHA)
    detection_control(nigoki, optimum)
    nigoki.save()

    with open(os.path.join(
        nigoki.outerpath, 'nigoki.pkl'
    ),'wb') as f: 
        pickle.dump(nigoki, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    )

    parser.add_argument(
        "-a", type=float, default=0.99, required=False,
        help="Alpha/significance parameter for detection control (default: 0.99)."
    )

    args = parser.parse_args()
    main(
        alpha=args.a,
    )