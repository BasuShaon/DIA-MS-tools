
#%%
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from sklearn.impute import KNNImputer

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

    data = data.dropna(axis=1, how='all')
    data = data.loc[:, ~(data == 0).all()]

    carrier.proteome = data
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
    means = data.mean(axis=0).replace(0, np.nan)

    carrier.pr_means_log2 = np.log2(means).replace([np.inf, -np.inf], np.nan)
    carrier.pr_miss_proportions = data.isna().mean(axis=0)

    return carrier


def detection_probability_curve(carrier, boundary=0.5):
    """
    Fits a logistic regression model for detection probability vs. log2 abundance.

    Parameters
    ----------
    carrier : object
        Must contain 'pr_means_log2' and 'pr_miss_proportions' attributes.
    boundary : float, optional
        Detection probability threshold for computing the decision boundary (default=0.5).

    Returns
    -------
    decision_boundary : float
        The abundance value at which detection probability equals the given boundary.
    result : GLMResults
        The fitted logistic regression model object.
    """
    df = pd.DataFrame({
        'Avg Log2 Abundance': carrier.pr_means_log2,
        'Detection Probability': 1 - carrier.pr_miss_proportions
    }).dropna()

    carrier.dpc_df = df.copy()

    X = sm.add_constant(df['Avg Log2 Abundance'])
    y = df['Detection Probability']

    model = sm.GLM(y, X, family=sm.families.Binomial())
    result = model.fit()

    intercept, coef = result.params
    decision_boundary = (np.log(boundary / (1 - boundary)) - intercept) / coef

    print(f"Decision boundary at {decision_boundary:.2f}")
    plot_detection_probability_curve(carrier, result, decision_boundary)

    return round(decision_boundary, 2), result


def plot_detection_probability_curve(carrier, result, decision_boundary):
    """
    Plots observed detection probabilities with the fitted logistic curve and decision boundary.

    Parameters
    ----------
    carrier : object
        Must contain dpc dataframe with 'Avg Log2 Abundance' and 'Detection Probability' columns.
    result : GLMResults
        The fitted logistic regression model.
    decision_boundary : float
        The abundance threshold where detection crosses the defined probability.

    Returns
    -------
    None. Displays a matplotlib plot.
    """
    df = carrier.dpc_df.copy()

    X_plot = np.linspace(df['Avg Log2 Abundance'].min(), df['Avg Log2 Abundance'].max(), 200)
    X_plot_const = sm.add_constant(X_plot)
    y_pred = result.predict(X_plot_const)

    plt.figure(plt.figure(figsize=(20, 10)))
    plt.scatter(df['Avg Log2 Abundance'], df['Detection Probability'], alpha=0.3, color='gray', label='Observed')
    plt.plot(X_plot, y_pred, linewidth=2, label='Logistic Fit')
    plt.axvline(decision_boundary, color='orange', linestyle='--', label='Decision Boundary')

    plt.xlabel('Avg Log2 Abundance')
    plt.ylabel('Detection Probability')
    plt.title('Detection Probability vs. Avg Log2 Abundance')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.savefig(os.path.join(carrier.outer_path, carrier.projectname + '_detectionprobabilitycurve.pdf'))
    plt.show()

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
        print(f"Dropping {fully_na_cols.size} fully NaN columns:", list(fully_na_cols))
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

  