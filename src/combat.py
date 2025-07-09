# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))

# third party imports

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc

def batch_correct(carrier):
    """
    Applies log2 transformation, ComBat batch correction to proteomics data and 
    visualizes coefficient of variation (CV) before and after correction using boxplots.

    Parameters:
    -----------
    carrier : object
        Object with `.proteome` DataFrame and `.metadata` DataFrame attributes.
        Assumes `metadata` contains 'MS.Batch' for batch assignment.
    toPlot : bool, optional (default=True)
        If True, generates CV boxplots before and after correction.

    Returns:
    --------
    outputdata : pd.DataFrame
        Batch-corrected proteomics data with original MultiIndex columns restored.
    """
    # Align metadata index
    carrier.reindex()

    # Print value counts by batch 
    print("Value counts for each batch:\n", carrier.metadata['MS.Batch'].value_counts().sort_index())

    batchdata = carrier.metadata['MS.Batch']
    inputdata = np.log2(carrier.proteome + 1)

    carrier.proteome_log2_beforecombat = inputdata.copy()

    # Save and flatten MultiIndex
    carrier.columns_multiindex = inputdata.columns
    inputdata.columns = ['|'.join(map(str, col)) for col in inputdata.columns]

    # Prepare AnnData object for ComBat
    obs = pd.DataFrame({'batch': batchdata}, index=inputdata.index)
    adata = sc.AnnData(
        X=inputdata.values,
        obs=obs,
        var=pd.DataFrame(index=inputdata.columns)
    )

    # Apply ComBat correction
    sc.pp.combat(adata, key='batch')

    # Convert corrected data back to DataFrame with original MultiIndex
    outputdata = pd.DataFrame(adata.X, index=inputdata.index, columns=carrier.columns_multiindex, dtype=float)
    
    carrier.proteome = outputdata
    carrier.status = carrier.status + '_combat'

    return carrier


def CV_plots(carrier, frame, title):
    """
    Computes and plots the coefficient of variation (CV) for each precursor grouped by batch.

    A boxplot of CVs across batches is generated and saved as a PDF.

    Parameters:
    -----------
    df : pd.DataFrame
        Proteomics data (e.g., log2-transformed intensities) with samples as rows and precursors as columns.
    
    batchdata : pd.Series
        Series of batch labels for each sample (index-aligned with `df`).
    
    title : str
        Title suffix used in the plot and filename.
    
    output_path : str
        Directory where the PDF plot will be saved.

    Returns:
    --------
    None
    """

    # Combine expression data with batch labels
    data = frame.copy()
    data['batchdata'] = carrier.metadata['MS.Batch'].copy()

    # Calculate CVs for each batch group
    def compute_cv(df):
        return df.std() / df.mean() * 100

    grouped = data.groupby('batchdata')

    cv_by_batch = grouped.apply(compute_cv)

    cv_by_batch.drop('batchdata', axis = 1, inplace = True)

    long_df = cv_by_batch.reset_index().melt(id_vars='batchdata', var_name='Precursor', value_name='CV')

    print(long_df)

    # Plot styling
    plt.rcParams.update({
        'axes.labelsize': 30,
        'axes.titlesize': 30,
        'legend.fontsize': 28,
        'xtick.labelsize': 26,
        'ytick.labelsize': 26
    })

    # Generate boxplot
    plt.figure(figsize=(20, 10))
    sns.boxplot(x='batchdata', y='CV', data=long_df, showfliers=False)
    plt.xticks(rotation=45)
    plt.xlabel('Batch')
    plt.ylabel('Coefficient of Variation (%)')
    plt.title(f'Distribution of CVs for Each Precursor Across Batches: {title}')
    plt.tight_layout()

    # Save plot
    plt.savefig(os.path.join(carrier.outer_path, f"{title.replace(' ', '_')}.pdf"))
    plt.show()