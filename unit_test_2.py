# %% Establish absolute path
import sys
import os
sys.path.append(os.path.abspath("src"))

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from carrier import Carrier
import detectioncontrol

# %%


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
