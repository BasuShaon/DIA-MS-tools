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
importlib.reload(combat)

d4 = combat.batch_correct(d3)
#d4.save()

# %%
importlib.reload(combat)
combat.CV_plots(d3, 'before Combat')
combat.CV_plots(d4, 'after Combat')

# %%

def check(df):
    print(f'min: {df.min().min()}')
    print(f'median: {df.median().median()}')
    print(f'zeroes: {(df == 0).sum().sum()}')
    print(f'na: {df.isna().sum().sum()}')
    print(f'inf: {np.isinf(df).sum().sum()}')

check(d3.proteome)

check(np.log2(d3.proteome))


# %%

sample = d3.proteome.sample(n = 6)

check(sample)

sample = pd.Series(sample.values.flatten())

sample.plot(kind = 'kde')
# %%

# Combine expression data with batch labels
data = d4.proteome.copy()
data['batchdata'] = d4.metadata['MS.Batch']

# %%

# Calculate CVs for each batch group
def compute_cv(df):
    return df.std() / df.mean() * 100

grouped = data.groupby('batchdata')

cv_by_batch = grouped.apply(compute_cv)

cv_by_batch.drop('batchdata', axis = 1, inplace = True)

long_df = cv_by_batch.reset_index().melt(id_vars='batchdata', var_name='Precursor', value_name='CV')


# %%
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

# %%
# Save plot
plt.savefig(os.path.join(output_path, f"{title.replace(' ', '_')}.pdf"))
plt.show()

# %%
print("long_df shape:", long_df.shape)
print(long_df.head())
print("Dtypes:\n", long_df.dtypes)
print("Nulls in batchdata:", long_df['batchdata'].isnull().sum())
print("Nulls in CV:", long_df['CV'].isnull().sum())
print("Unique batches:", long_df['batchdata'].unique())
print("All CVs numeric:", pd.to_numeric(long_df['CV'], errors='coerce').notnull().all())
