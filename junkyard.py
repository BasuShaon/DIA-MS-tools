# %% Establish absolute path
import sys
import os
sys.path.append(os.path.abspath("src"))
from carrier import Carrier
import detectioncontrol
import imputer
import combat
import precursor2protein
import pandas as pd

# %% Clean up and format precursor and metadata 

data_dir = os.path.join(os.path.dirname(__file__), 'input')
print(data_dir)

secondpass = pd.read_csv(os.path.join(data_dir, 'SB_proteome-druggable-genome_KO_prmatrix-secondpass_250320a.csv'), index_col = [1,2]).iloc[:, 1:]
secondpass.columns = secondpass.columns.str.split('/').str[-1]

metadata = pd.read_csv(os.path.join(data_dir, 'SB_proteome-druggable-genome_KO_metadata_251007a.csv'), 
        index_col = 1, sep = ',').dropna(axis = 1, how = 'all')

# %%

secondpass_dr = pd.read_csv(
    os.path.join(
        data_dir, 
        'SB_proteome-druggable-genome_DR_prmatrix-secondpass_250320a.csv'), index_col = [1,2]).iloc[:, 1:]


metadata_dr = pd.read_csv(os.path.join(data_dir, 'SB_proteome-druggable-genome_DR_metadata_251007a.csv'), 
        index_col = 0, sep = ',').dropna(axis = 1, how = 'all')

# %% Create Carrier object 
shogoki = Carrier(proteome=secondpass, metadata=metadata, outerpath='/Users/shaon/Desktop/50_0102/data', projectname='SB_500102_secondpass')                

# %% Detection Control 
optimum = detectioncontrol.calculate_optimum_threshold(shogoki, alpha = 0.99)
detectioncontrol.detection_control(shogoki, optimum)
shogoki.save()

# %% Mixed Imputation 
imputer.preprocess_data(shogoki)
imputer.compute_log2_means_and_missingness(shogoki)
bound, glm = imputer.detection_probability_curve(shogoki)
imputer.mixed_imputation_in_batch(shogoki, bound)
shogoki.save()

# %% Combat Correction
combat.batch_correct(shogoki)
combat.CV_plots(shogoki, shogoki.proteome_log2_beforecombat, 'before Combat')
combat.CV_plots(shogoki, shogoki.proteome, 'after Combat')
shogoki.save()

# %% MaxLFQ 
precursor2protein.format_df_maxlfq(shogoki)
shogoki.save()
precursor2protein.maxlfq(shogoki, 
                         os.path.join(shogoki.outer_path, 
                                      shogoki.projectname + '_' + 
                                      shogoki.status + '_' + 
                                      shogoki.thedate + '.tsv'))
shogoki.status

# %% Detection Control Drop out Analysis -> Add into Dcontrol Module

#boundary = 82.12
#miss = secondtype.missingness.copy()
#miss.index = miss.index.str.split('/').str[-1]
#miss.head()
#miss.loc[miss<=boundary].describe()
#miss.loc[miss>boundary].describe()

#miss = miss.to_frame()  
#miss['Category'] = np.where(miss.iloc[:,0] <= boundary, 'Retain', 'Drop')
#plt.figure(figsize=[3,6])
#sns.violinplot(data=miss, y=0, hue='Category', fill=True, alpha=0.4)
#plt.title('KDE of Missing Distributions')
#plt.xlabel('Value')
#plt.axhline(boundary, linestyle = '--', color = 'red')
#plt.ylabel('Density')
#plt.savefig('/Users/shaon/Desktop/50_0102/data/SB_500102_secondpass_dcontrol_82_splitkde_250513.png')
#plt.show()
# 
#miss = secondtype.missingness.copy()
#miss.index = miss.index.str.split('/').str[-1]

#batch_comments = batch['Comment'].copy()
#batch_comments.fillna('no technical', inplace = True)

#subset_comments = batch_comments.loc[miss>=82.12]
#batch_comments = batch_comments.loc[miss<82.12]

#all_categories = sorted(set(batch_comments.unique()) | set(subset_comments.unique()))
#palette = sns.color_palette('hsv_r', n_colors=len(all_categories))
#color_map = dict(zip(all_categories, palette))
#batch_counts = batch_comments.value_counts()
#batch_colors = [color_map[label] for label in batch_counts.index]

#batch_counts.plot(
#    kind='pie',
#    autopct='%1.1f%%',
#    startangle=90,
#    figsize=(10, 10),
#    textprops={'fontsize': 16},
#    colors=batch_colors
#)
#plt.ylabel('')
#plt.tight_layout()
#plt.show()
#subset_counts = subset_comments.value_counts()
#subset_colors = [color_map[label] for label in subset_counts.index]

#subset_counts.plot(
#    kind='pie',
#    autopct='%1.1f%%',
#    startangle=90,
#    figsize=(10, 10),
#    textprops={'fontsize': 16},
#    colors=subset_colors
#)
#plt.ylabel('')
#plt.tight_layout()
#plt.show()
