# %% Establish absolute path
import sys
import os
sys.path.append(os.path.abspath("src"))
import importlib
import pandas as pd
import numpy as np 
import seaborn as sns
import matplotlib.pyplot as plt
from carrier import Carrier
import detectioncontrol
import imputer
import combat
import precursor2protein

# %% Detection Control

metapath = '/Users/shaon/Desktop/50_0102/2025-02-03_AF_50-0102_metadata_with_DIANN_output.txt'

firstpass = pd.read_csv('/Users/shaon/Desktop/50_0102/data/SB_500102_firstpass_250320a.tsv', index_col = [1,2]).iloc[:, 1:]

secondpass = pd.read_csv('/Users/shaon/Desktop/50_0102/data/SB_500102_secondpass_250320a.tsv', index_col = [1,2]).iloc[:, 1:]

batch = pd.read_csv('/Users/shaon/Desktop/50_0102/2025-02-03_AF_50-0102_metadata_with_DIANN_output.txt', 
        index_col = 0, sep = '\t').dropna(axis = 1, how = 'all')


zerotype = Carrier(proteome=firstpass, metadata=batch['MS.Batch'], outerpath='/Users/shaon/Desktop/50_0102/data', projectname='SB_500102_firstpass')                

#secondtype= Carrier(proteome=secondpass, metadata=batch['MS.Batch'], outerpath='/Users/shaon/Desktop/50_0102/data', projectname='SB_500102_secondpass')                

def dcontrol(inputxx):
    opt99 = detectioncontrol.calculate_optimum_threshold(inputxx, alpha = 0.99)
    #opt95 = detectioncontrol.calculate_optimum_threshold(inputxx, alpha = 0.95)
    detectioncontrol.detection_control(inputxx,opt99)
    inputxx.save()

dcontrol(zerotype)

# %%

miss = zerotype.missingness.copy()
miss.index = miss.index.str.split('/').str[-1]
miss.head()

#%%
miss.loc[miss<=85.97].describe()
# %%
miss.loc[miss>85.97].describe()

# %%

batch_comments = batch['Comment'].copy()
batch_comments.fillna('no technical', inplace = True)

subset_comments = batch_comments.loc[miss>=85.97]
batch_comments = batch_comments.loc[miss<85.97]

# %%
all_categories = sorted(set(batch_comments.unique()) | set(subset_comments.unique()))
palette = sns.color_palette('hsv_r', n_colors=len(all_categories))
color_map = dict(zip(all_categories, palette))
batch_counts = batch_comments.value_counts()
batch_colors = [color_map[label] for label in batch_counts.index]

batch_counts.plot(
    kind='pie',
    autopct='%1.1f%%',
    startangle=90,
    figsize=(10, 10),
    textprops={'fontsize': 16},
    colors=batch_colors
)
plt.ylabel('')
plt.tight_layout()
plt.show()
subset_counts = subset_comments.value_counts()
subset_colors = [color_map[label] for label in subset_counts.index]

subset_counts.plot(
    kind='pie',
    autopct='%1.1f%%',
    startangle=90,
    figsize=(10, 10),
    textprops={'fontsize': 16},
    colors=subset_colors
)
plt.ylabel('')
plt.tight_layout()
plt.show()


# %% Impute

data = pd.read_csv('/Users/shaon/Desktop/50_0102/data/SB_500102_secondpass_dcontrol_82_250513.tsv', index_col = [0,1]).T

data.index = data.index.str.replace('/s-mcpb-ms03/union/is/PRO1/Data/50-0102/', '', regex=False)

metadata = pd.read_csv('/Users/shaon/Desktop/50_0102/2025-02-03_AF_50-0102_metadata_with_DIANN_output.txt', 
        index_col = 0, sep = '\t').dropna(axis = 1, how = 'all')

metadata = metadata.loc[data.index]

test_type = Carrier(data, metadata, '/Users/shaon/Desktop/50_0102/data', 'SB_500102_secondpass')                

# %%

importlib.reload(imputer)
importlib.reload(combat)

d1 = imputer.preprocess_data(test_type)
d2 = imputer.compute_log2_means_and_missingness(d1)
bound, result = imputer.detection_probability_curve(d2)

d3 = imputer.mixed_imputation_in_batch(d2, bound)

#d3.save()


# %%

d4 = combat.batch_correct(d3)

combat.CV_plots(d4, d4.proteome_log2_beforecombat, 'before Combat')

combat.CV_plots(d4, d4.proteome, 'after Combat')

#d4.save()

# %%

def check(df):
    print(f'min: {df.min().min()}')
    print(f'median: {df.median().median()}')
    print(f'zeroes: {(df == 0).sum().sum()}')
    print(f'na: {df.isna().sum().sum()}')
    print(f'inf: {np.isinf(df).sum().sum()}')


# %% Summarization module

data = pd.read_csv('/Users/shaon/Desktop/50_0102/data/SB_500102_secondpass_dcontrol_imputed_combat_250704.tsv', 
                   header = [0,1], index_col = [0])


# %%
importlib.reload(precursor2protein)

d5 = precursor2protein.format_df_maxlfq(d4)

d5.save()

# %%
importlib.reload(precursor2protein)
# Update long form filepath into module

precursor2protein.maxlfq(d5, '/Users/shaon/Desktop/50_0102/data/SB_500102_secondpass__imputed_inbatch_combat_long_250707.tsv')
