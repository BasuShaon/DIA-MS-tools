from datetime import datetime
import pandas as pd
import numpy as np
import os
import pickle
import argparse

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
        self.outerpath = outerpath
        self.projectname = projectname

        #detection control attrs
        self.missingness = None

        #mixed imputation attrs
        self.pr_means_log2 = None
        self.pr_miss_proportions = None

        #batch correction attrs
        self.proteome_log2_beforecombat = None
        self.columns_multiindex = None

        #maxlfq r comms attrs
        self.thedate = datetime.now().strftime("%y%m%d")
        

    def save(self):
        """  
        Saves proteome as a .tsv with a specific, dated file name.

        Parameters
        ----------      
            filename : str
                base prefix, i.e. '/Desktop/SB_project1_prmatrix'

        """

        # construct savepath
        savename = self.projectname + '_' + self.status + '_' + self.thedate + '.tsv'
        print(savename)

        self.proteome.to_csv(os.path.join(self.outerpath, savename))

        return savename
    
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

def init_carrier(data_df, meta_df, proj):
    # Configuration
    DATA_DF = data_df
    META_DF = meta_df
    PROJECT_NAME = proj
    DATA_DIR = os.path.join(os.path.dirname(__file__), '..', 'input')
    OUT_DIR = os.path.join(os.path.dirname(__file__), '..', 'output')

    print("Starting proteomics pipeline...")
    
    print("\n1. Loading data...")
    secondpass = pd.read_csv(
        f'{DATA_DIR}/{DATA_DF}',
        index_col=[1, 2]
    ).iloc[:, 1:]
    secondpass.columns = secondpass.columns.str.split('/').str[-1]
    
    metadata = pd.read_csv(
        f'{DATA_DIR}/{META_DF}',
        index_col=1,
        sep=','
    ).dropna(axis=1, how='all')
    
    print("2. Initializing Carrier...")
    shogoki = Carrier(
        proteome=secondpass,
        metadata=metadata,
        outerpath=OUT_DIR,
        projectname=PROJECT_NAME
    )
    
    with open(os.path.join(
        shogoki.outerpath, 'shogoki.pkl'
    ),'wb') as f: 
        pickle.dump(shogoki, f)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
    )

    parser.add_argument(
        "-d", required=True,
        help="Filename of the raw data CSV (relative to './input')."
    )
    parser.add_argument(
        "-m", required=True,
        help="Filename of the metadata CSV (relative to './input')."
    )
    parser.add_argument(
        "-p", required=True,
        help="Project name used as output prefix."
    )

    args = parser.parse_args()
    init_carrier(
        data_df=args.d,
        meta_df=args.m,
        proj=args.p,
    )