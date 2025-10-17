
#%%
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.impute import KNNImputer
import argparse
import pickle
from carrier import Carrier

def preprocess_data(carrier):
    """
    Removes features (columns) from the proteome data that contain only NaNs or only zeros.

    Parameters
    ----------
    carrier : object
        An object with a 'proteome' attribute (pd.DataFrame) containing raw intensity data.

    Returns
    -------
    carrier : object
        The input object with its 'proteome' attribute cleaned.
    """
    data = carrier.proteome.copy()

    n_all_nan = data.isna().all().sum()
    n_all_zero = (data == 0).all().sum()
    print(f"Removing {n_all_nan} features with all NaNs")
    print(f"Removing {n_all_zero} features with all zeros")

    data = data.dropna(axis=1, how='all') # drop all-NaN columns
    data = data.loc[:, ~(data == 0).all()] # drop all-zero columns

    carrier.proteome = data.T # transpose proteome in carrier object for all downstream modules 
    return carrier

def compute_log2_means_and_missingness(carrier):
    """
    Calculates per-feature log2 mean intensities and missingness proportions.

    Parameters
    ----------
    carrier : object
        An object with a 'proteome' attribute (pd.DataFrame).

    Returns
    -------
    carrier : object
        The same object with new attributes:
            - pr_means_log2: pd.Series of log2-transformed mean values.
            - pr_miss_proportions: pd.Series of missing value proportions.
    """
    data = carrier.proteome.copy() 

    means = data.mean(axis=0).replace(0, np.nan) # replace zero means (by sample) to NaN beofre the log2 transformation

    carrier.pr_means_log2 = np.log2(means).replace([np.inf, -np.inf], np.nan) # means log2, with infinities removed
    carrier.pr_miss_proportions = data.isna().mean(axis=0) # average detection proportion by sample
    carrier.pr_miss_counts = data.isna().sum(axis=0) # detection counts by samples

    return carrier

def detection_probability_curve(carrier, boundary=0.5):
    """
    Fits logistic regression for detection probability vs. log2 abundance,
    using response-binned data (100 bins).

    Parameters
    ----------
    carrier : object
        Must contain 'pr_means_log2' and 'pr_miss_proportions'.
    boundary : float, optional
        Detection probability threshold for decision boundary (default=0.5).

    Returns
    -------
    decision_boundary : float
        Abundance at the given detection probability boundary.
    result : GLMResults
        Logistic regression model fit on binned response data.

    Author SB 
    
    """
    # Raw data preparation
    df = pd.DataFrame({
        'Avg Log2 Abundance': carrier.pr_means_log2,
        'Detection Probability': 1 - carrier.pr_miss_proportions
    }).dropna()

    # Bin the response variable (Detection Probability)
    df['response_bin'] = pd.qcut(df['Detection Probability'], q=100, duplicates='drop')

    # Compute mean abundance and mean detection probability per bin
    grouped = df.groupby('response_bin').agg(
        abundance_mean=('Avg Log2 Abundance', 'median'),
        detection_mean=('Detection Probability', 'median')
    ).dropna()

    # Prepare data for GLM
    X = sm.add_constant(grouped['abundance_mean'])
    y = grouped['detection_mean']

    # Fit binomial GLM (using proportions directly with freq_weights)
    model = sm.GLM(y, X, family=sm.families.Binomial())
    result = model.fit()

    # Decision boundary calculation
    intercept, coef = result.params
    decision_boundary = (np.log(boundary / (1 - boundary)) - intercept) / coef

    # Save summary dataframe to carrier
    carrier.dpc_df = grouped.rename(columns={
        'abundance_mean': 'Avg Log2 Abundance',
        'detection_mean': 'Detection Probability'
    })

    print(f"Decision boundary at {decision_boundary:.2f}")
    plot_detection_probability_curve(carrier, result, decision_boundary)

    return round(decision_boundary, 2), result

def plot_detection_probability_curve(carrier, result, decision_boundary):
    """
    Plots logistic regression fit (from response-binned data) over original raw data.

    Parameters
    ----------
    carrier : object
        Must contain 'pr_means_log2', 'pr_miss_proportions', and 'dpc_df'.
    result : GLMResults
        Logistic regression fit on binned data.
    decision_boundary : float
        Abundance threshold at specified detection probability.

    Returns
    -------
    None. Displays plot and saves PDF.
    """
    # Raw data
    df_raw = pd.DataFrame({
        'Avg Log2 Abundance': carrier.pr_means_log2,
        'Detection Probability': 1 - carrier.pr_miss_proportions
    }).dropna()

    # Generate smooth logistic curve for plotting
    X_plot = np.linspace(df_raw['Avg Log2 Abundance'].min(), df_raw['Avg Log2 Abundance'].max(), 300)
    X_plot_const = sm.add_constant(X_plot)
    y_pred = result.predict(X_plot_const)

    # Plotting setup
    plt.figure(figsize=(12, 7))

    # Scatter raw data
    plt.scatter(
        df_raw['Avg Log2 Abundance'], 
        df_raw['Detection Probability'], 
        alpha=0.1, color='gray', s=10, label='Raw Data'
    )

    # Binned averages (for clarity)
    #plt.scatter(
    #    carrier.dpc_df['Avg Log2 Abundance'], 
    #    carrier.dpc_df['Detection Probability'], 
    #    color='blue', s=50, alpha=0.8, label='Binned Means'
    #)

    # Logistic fit line from binned response fit
    plt.plot(
        X_plot, y_pred, color='blue', linewidth=2.5, label='Logistic Fit (Binned Response)'
    )

    # Decision boundary line
    plt.axvline(
        decision_boundary, color='orange', linestyle='--', linewidth=2,
        label=f'Decision Boundary @ {decision_boundary:.2f}'
    )

    # Finalize plot
    plt.xlabel('Avg Log2 Abundance', fontsize=14)
    plt.ylabel('Detection Probability', fontsize=14)
    plt.title('Detection Probability Curve (Logistic Fit from Response-Binned Data)', fontsize=16)
    plt.legend(fontsize=12)
    plt.grid(True, linestyle='--', alpha=0.5)
    plt.tight_layout()

    # Save plot
    plt.savefig(
        os.path.join(carrier.outerpath, carrier.projectname + '_detectionprobabilitycurve.pdf'),
        dpi=300
    )
    # plt.show()


def mixed_imputation_global(carrier, boundary, knn=3):
    """
    Performs mixed imputation on proteomics data, using:
      - Global minimum imputation for MNAR features
      - Global KNN imputation for MAR features (no batch separation)

    Parameters:
    -----------
    carrier : object
        An object with the following required attributes:
            - carrier.proteome : pd.DataFrame
                Raw proteomics data with samples as rows and features (precursors) as columns.
            - carrier.metadata : pd.DataFrame
                Metadata including a 'MS.Batch' column indicating batch membership.
            - carrier.pr_means_log2 : pd.Series
                Log2-transformed global mean intensities of features; used to classify MNAR/MAR.
    boundary : float
        Threshold used to classify features as MNAR (if mean < boundary) or MAR (if >= boundary).
    knn : int, optional (default=3)
        Number of neighbors to use in global KNN imputation for MAR features.

    Returns:
    --------
    carrier : object
        The input object with updated:
            - carrier.proteome : pd.DataFrame
                Imputed data with MNAR and MAR features appropriately filled.
            - carrier.status : str
                Set to 'dcontrol_imputed' to indicate imputation is complete.

    Notes:
    ------
    - MNAR features are imputed using global minimum value imputation.
    - MAR features are imputed using global KNN imputation (ignoring batch information).
    - Columns fully missing (NaN) are dropped.
    """

    data = carrier.proteome.copy()

    # Identify MNAR and MAR columns
    mnar_cols = carrier.pr_means_log2[carrier.pr_means_log2 < boundary].index.intersection(data.columns)
    mar_cols = data.columns.difference(mnar_cols)

    # MNAR imputation with global minimum
    min_values = data[mnar_cols].min()
    data[mnar_cols] = data[mnar_cols].fillna(min_values)

    # Check and drop fully NaN MAR columns
    fully_na_cols = data[mar_cols].columns[data[mar_cols].isna().all()]
    if len(fully_na_cols) > 0:
        print(f"\nDropping {fully_na_cols.size} fully NaN columns:", list(fully_na_cols))
        data = data.drop(columns=fully_na_cols)
        mar_cols = mar_cols.difference(fully_na_cols)

    # Global KNN imputation for MAR features
    imputer = KNNImputer(n_neighbors=knn)
    mar_imputed_array = imputer.fit_transform(data[mar_cols])

    mar_imputed_df = pd.DataFrame(
        mar_imputed_array,
        index=data.index,
        columns=mar_cols
    )

    # Re-combine MNAR and MAR imputed columns
    data_imputed = pd.concat([data[mnar_cols], mar_imputed_df], axis=1)[data.columns]

    carrier.proteome = data_imputed
    carrier.status = carrier.status + '_imputed_global'

    return carrier

def mixed_imputation_in_batch(carrier, boundary, knn=3):
    """
    Performs mixed imputation on proteomics data, using:
      - Global minimum imputation for MNAR features
      - Per-batch KNN imputation for MAR features

    Parameters:
    -----------
    carrier : object
        An object with the following required attributes:
            - carrier.proteome : pd.DataFrame
                Raw proteomics data with samples as rows and features (precursors) as columns.
            - carrier.metadata : pd.DataFrame
                Metadata including a 'MS.Batch' column indicating batch membership.
            - carrier.pr_means_log2 : pd.Series
                Log2-transformed global mean intensities of features; used to classify MNAR/MAR.
    boundary : float
        Threshold used to classify features as MNAR (if mean < boundary) or MAR (if >= boundary).
    knn : int, optional (default=3)
        Number of neighbors to use in KNN imputation for MAR features.

    Returns:
    --------
    carrier : object
        The input object with updated:
            - carrier.proteome : pd.DataFrame
                Imputed data with MNAR and MAR features appropriately filled.
            - carrier.status : str
                Set to 'dcontrol_imputed' to indicate imputation is complete.

    Notes:
    ------
    - MNAR features are imputed using global minimum value imputation.
    - MAR features are imputed using KNN within each batch.
    - Any columns that are fully missing (NaN) in a given batch are dropped across all batches.
    - The function assumes that batch effects may influence missingness patterns and treats each batch separately.
    """
    
    data = carrier.proteome.copy()
    metadata = carrier.metadata.copy()

    #re-index
    metadata = metadata.loc[data.index]

    # Identify MNAR columns based on the boundary
    mnar_cols = carrier.pr_means_log2[carrier.pr_means_log2 < boundary].index.intersection(data.columns)
    mar_cols = data.columns.difference(mnar_cols)

    # Impute MNAR features globally (minimum imputation)
    min_values = data[mnar_cols].min()
    data[mnar_cols] = data[mnar_cols].fillna(min_values)

    # KNN imputation within each batch for MAR features
    imputed_batches = []
    cols_to_drop = set()

    for batch_id in metadata['MS.Batch'].unique():
        batch_samples = metadata[metadata['MS.Batch'] == batch_id].index
        batch_data = data.loc[batch_samples]

        # Select only MAR columns for KNN
        mar_batch_data = batch_data[mar_cols]

        # Identify fully NaN columns in this batch
        fully_na_cols = mar_batch_data.columns[mar_batch_data.isna().all()]

        if len(fully_na_cols) > 0:
            print(f"Batch '{batch_id}' - Dropping {fully_na_cols.size} fully NaN columns:", list(fully_na_cols))
            cols_to_drop.update(fully_na_cols)

        # Drop columns fully NaN within this batch
        mar_batch_data_non_na = mar_batch_data.drop(columns=fully_na_cols)

        # Perform KNN imputation
        imputer = KNNImputer(n_neighbors=knn)
        mar_imputed_array = imputer.fit_transform(mar_batch_data_non_na)

        mar_imputed_df = pd.DataFrame(
            mar_imputed_array,
            index=batch_samples,
            columns=mar_batch_data_non_na.columns
        )

        # Re-combine MNAR and MAR imputed columns (excluding fully NaN dropped columns)
        batch_imputed = pd.concat(
            [batch_data[mnar_cols], mar_imputed_df],
            axis=1
        )

        imputed_batches.append(batch_imputed)

    # Drop columns fully NaN in any batch from all batches
    if cols_to_drop:
        print(f"Dropping {len(cols_to_drop)} columns fully NaN in at least one batch:", list(cols_to_drop))
        imputed_batches = [df.drop(columns=cols_to_drop, errors='ignore') for df in imputed_batches]

    # Combine imputed batches
    data_imputed = pd.concat(imputed_batches).loc[data.index]

    carrier.proteome = data_imputed
    carrier.status = carrier.status + '_imputed_inbatch'

    return carrier


def main(bound):
    # Configuration
    BOUND = bound

    with open(os.path.join(os.path.dirname(__file__), "../output", "nigoki.pkl"), "rb") as f:
        sangoki = pickle.load(f)

    print("4. Running imputation...")
    preprocess_data(sangoki)
    compute_log2_means_and_missingness(sangoki)
    if bound == 'Auto':
        bound, glm = detection_probability_curve(sangoki)
        mixed_imputation_in_batch(sangoki, bound)
    else:
        bound = float(bound)
        mixed_imputation_in_batch(sangoki, bound)
    sangoki.save()

    with open(os.path.join(
        sangoki.outerpath, 'sangoki.pkl'
    ),'wb') as f: 
        pickle.dump(sangoki, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    )

    parser.add_argument(
        "-b", default="Auto",
        help="Boundary for MNAR/MAR split. Float or 'Auto' to fit via logistic curve (default: Auto)."
    )

    args = parser.parse_args()
    main(
        bound=args.b,
    )
# %%
