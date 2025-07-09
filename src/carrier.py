from datetime import datetime
import pandas as pd
import numpy as np
import os

class Carrier: 
    
    def __init__(self, proteome, metadata, outerpath, projectname):
        """  
        Parameters
        ----------      
            proteome : pd.DataFramee
                Precursor intensity data
            metadata : pd.Series
                Associated metadata
            outerpath : str
                outpath to save files
            projectname : str
                common project suffix
        """

        self.status = ''
        self.proteome = proteome
        self.metadata = metadata
        self.outer_path = outerpath
        self.projectname = projectname

        #detection control attrs
        self.missingness = None

        #mixed imputation attrs
        self.pr_means_log2 = None
        self.pr_miss_proportions = None

        #batch correction attrs
        self.proteome_log2_beforecombat = None
        self.columns_multiindex = None

        #maxlfq attrs

        

    def save(self):
        """  
        Saves proteome as a .tsv with a specific, dated file name.

        Parameters
        ----------      
            filename : str
                base prefix, i.e. '/Desktop/SB_project1_prmatrix'
        """

        # construct savepath
        thedate = datetime.now().strftime("%y%m%d")
        savename = self.projectname + '_' + self.status + '_' + thedate + '.tsv'
        print(savename)

        self.proteome.to_csv(os.path.join(self.outer_path, savename))

        return 
    
    def reindex(self):
        """ 
        Slice metadata by proteome index.
        """  

        self.metadata = self.metadata.loc[self.proteome.index]

    def describe(self):
        """ 
        Quick df_describe() on the current carried proteome df
        """  
        df = self.proteome.copy()

        print(f'min: {df.min().min()}')
        print(f'median: {df.median().median()}')
        print(f'zeroes: {(df == 0).sum().sum()}')
        print(f'na: {df.isna().sum().sum()}')
        print(f'inf: {np.isinf(df).sum().sum()}')

