from datetime import datetime
import pandas as pd
import numpy as np
import os

class Carrier: 
    """ 
    Data class that contains the proteome and related processing information. 

    Attributes
    ----------
        status : str
            status of the carrier frame
        proteome : pd.DataFrame
            precursor intensity data 
        batch : pd.Series
            batch data
        save(filename) : function
            saves frame in its current state 
    """
    
    def __init__(self, proteome, metadata, outerpath, projectname):
        """  
        Parameters
        ----------      
            intensity_data : pd.DataFramee
                Precursor intensity data
            batch_data : pd.Series
                MS Batch data 
        """
        self.status = None
        self.proteome = proteome
        self.metadata = metadata
        self.outer_path = outerpath
        self.projectname = projectname

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
    
    def save_joblib():
        return
