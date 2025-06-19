# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import ListedColormap
import seaborn as sns
from sklearn.impute import KNNImputer
import statsmodels.api as sm


# %%

def detectionprobabilitycurve(data_in, boundary, ramp):

        
    #drop missing in all columns
    data_in.dropna(axis = 1, how = 'all', inplace = True)

    print('Replace ' + str((data_in == 0).sum().sum()) + '  features with 100% zero abundance values / sample with NaN')


    # NaN to 0 
    data_in.replace(0, np.nan, inplace = True)

    # L2 transform

      # calculate precursor stats for DPC curve

    means = data_in.mean(axis=0)

    means_log2 = np.log2(means.astype(float))

    means_log2.replace([np.inf, -np.inf], np.nan, inplace=True)

    missing = data_in.isnull().mean(axis=0)

    missing_df = pd.DataFrame({
        'Avg Precursor Abundance' : means_log2.values, 
        'Missing %' : 1-(missing.values)})

    # Remove NA, Fit model

    cleaned_df = missing_df.dropna(subset=['Avg Precursor Abundance'])

    X = cleaned_df[['Avg Precursor Abundance']]

    y = cleaned_df['Missing %']
            
    # Add a constant to the model for the intercept
    X = sm.add_constant(X)

    # Fit the Logistic Regression Model using GLM
    model = sm.GLM(y, X, family=sm.families.Binomial(sm.families.links.logit()))

    result = model.fit()

    # Extract the intercept and coefficient
    intercept, coef = result.params

    result = model.fit()
    # Calculate expected value of missing proportions
    # Remove NA
    # Fit the GLM 
    # Plot and extract boundary and stats of the model 

    avg_precursor_abundance = (np.log(boundary / (1 - boundary)) - intercept) / coef

    boundary = round_up(avg_precursor_abundance, 2) 

    print(result.summary())

    mnar = cleaned_df[cleaned_df['Avg Precursor Abundance'] < boundary]

    mar = cleaned_df[cleaned_df['Avg Precursor Abundance'] >= boundary]

    omit = cleaned_df[cleaned_df['Missing %'] <= ramp/100]

    # Generating a range of values for Avg Precursor Abundance

    x_range = np.linspace(X['Avg Precursor Abundance'].min(), X['Avg Precursor Abundance'].max(), 100)

    x_range = sm.add_constant(pd.DataFrame(x_range, columns=['Avg Precursor Abundance']))

    # Predicting probabilities

    y_pred = result.predict(x_range)

    plt.figure(figsize=(12, 10))

    plt.rcParams['axes.labelsize'] = '30'   

    plt.rcParams['axes.titlesize'] = '30' 

    plt.rcParams['legend.fontsize'] = '28'  

    plt.rcParams['xtick.labelsize'] = '26'  

    plt.rcParams['ytick.labelsize'] = '26' 
    
    plt.plot(x_range['Avg Precursor Abundance'], y_pred, marker = 'o', label = 'p-val: ' + str(p))  

    plt.axvline(boundary, 
        label = 'Decision boundary (' + str(boundary*100) + '%): ' + str(boundary), 
        color = 'red', linestyle = '--')
    
    plt.axhline(ramp/100,
        color = 'black', linestyle = '--')

    plt.scatter(mar['Avg Precursor Abundance'],
                
        mar['Missing %'], alpha=0.1, s = 1, color = 'orange', label = 'knn(3)-impute')
        
    plt.scatter(mnar['Avg Precursor Abundance'],
            mnar['Missing %'], alpha=0.1, s = 1, color = 'yellow', label = 'min-impute')
    
    plt.scatter(omit['Avg Precursor Abundance'],
            omit['Missing %'], alpha=0.1, s = 1, color = 'grey', label = 'Below threshold (' + str(ramp) + '%): drop')
    
    # Labeling the plot

    plt.xlabel('Average Log2 Precursor Abundance')

    plt.ylabel('Average Precursor Decection')

    plt.title('Detection Probability Curve')   

    print('Decision boundary value (avg log2)= ' + str(boundary))

    plt.grid(True)

    plt.legend()

    plt.show()

    return 
