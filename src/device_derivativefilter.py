# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
from device import Device

# third party imports
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

class DerivativeFilter(Device):
    """ 
    Performs thresholding by modelling precursor detection across matrix dimensions. 
    """
    def __init__(self, data):
        """ 
        Parameters:
        -----------
        prmatrix : pd.DataFrame
            The precusor matrix (indexed by precursor, features samples)

        step : int
            number of intervals for thresholding between 0 - 100%  completeness,
            default = 20 

        retention_stats: pd.DataFrame
            df containing retention statistics 
        
        thresh_prmatrix: pd.DataFrame
            The thresholded prmatrix after retention max fitting

        """
        self.prmatrix = data
        self.step = None
        self.retention_stats = None
        self.maximum = None
        self.filtered_prmatrix = None
        return
    
    def apply_retention_filter(self, manual = None):
        """ 
        This main function applies the calculate_sample_stats(), then 
        find_optimal_thresholds(). 

        It logs the whole process during script runtime. 

        Then it uses the threshold to filter the pr_matrix, and returns it.

        Author: SB

        Returns
        -------
        filtered_prmatrix : pd.DataFrame
            Filtered DataFrame based on sample completeness if optimum is found, 
            if no minimum exists, returns original DataFrame without filtering on
            completeness 

        """
        df = self.prmatrix.copy()

        if manual is None:
            optimal_threshold = self.maximum
            # Plot retention curve
            sns.lineplot(x=retention_stats.index * (100 / step), 
                y=retention_stats['Missing Sample Percentage'])
            plt.axvline(x=optimal_threshold, linestyle="--", color="red", 
                label=f"Threshold: {optimal_threshold:.2f}%")
            plt.ylabel('Missing Sample Percentage')
            plt.xlabel('Sample Completeness (%)')
            plt.title('Retention Threshold Detection')
            plt.legend()
            plt.show()

        else: 
            optimal_threshold = manual

        # Apply filtering
        filtered_df = df[df.isnull().mean(axis=1) * 100 <= optimal_threshold]

        self.filtered_prmatrix  = filtered_df.copy()
        return self.filtered_prmatrix
    
    def calculate_retention_stats(self, step = 20):
        """ 
        Calculates number of recovered peptides across incremental threshold intervals 
        between 0 - 100% sample completeness.
        Author: SB
        
        step : int
            number of intervals for thresholding bteween 0 - 100%  completeness,
            default = 20 

        Returns:
        --------
        retention_stats : pd.DataFrame
            A DataFrame containing the retention statistics with:
            - Samples: Number of retained samples at each threshold.
            - Precursors: Number of retained protein features.
            - Missing Sample Percentage: Average percentage of missing values.

        """
        summary = []
        # Remove rows that are entierely NaN
        data = self.prmatrix.dropna(axis = 0, how = 'all')

        # Calculate retention thresholds
        for x in range(step + 1):  
            gradient = x * (100 / step)
            select_sample = data.isnull().mean(axis=1) * 100
            sample_filtered = data.loc[select_sample <= gradient, :]
            if sample_filtered.empty:
                sample_miss = 0.0
            else:
                sample_miss = sample_filtered.isnull().mean(axis=0).mean() * 100
            summary.append((sample_filtered.shape[0], sample_filtered.shape[1], sample_miss))

        # Return retention data
        self.retention_stats = pd.DataFrame(summary, columns=['Samples', 'Precursors', 
        'Missing Sample Percentage'])
        self.step = step

        return self.retention_stats
        
    def find_optimal_threshold(self):
        """ 
        Identifies the threshold at which missing sample percentage stabilizes 
        by finding first and second derivative minima on retention data

        Returns
        -------
        float or None
            Optimal completeness threshold or None if no stabilization point is found.

        """
        # Set threshold multiplier and calculate gradient minimum
        xmult = 100 / self.step
        samples = self.retention_stats['Missing Sample Percentage']
        threshold = self.retention_stats.index
        first_derivative = np.gradient(samples, threshold)
        second_derivative = np.gradient(first_derivative, threshold)
        zero_idx = np.where(np.isclose(second_derivative, 0, atol = 1e-6))[0][-1] #last zero
        if zero_idx.size > 0:
            self.maximum = zero_idx * mult


