from datetime import datetime
import pandas as pd
import numpy as np

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
    
    def __init__(self, intensity_data, batch_data, path):
        """  
        Parameters
        ----------      
            intensity_data : pd.DataFramee
                Precursor intensity data
            batch_data : pd.Series
                MS Batch data 
        """
        self.status = 'raw'
        self.proteome = intensity_data
        self.batch = batch_data
        self.outer_path = path

    def save(self, filename):
        """  
        Saves proteome as a .tsv with a specific, dated file name.

        Parameters
        ----------      
            filename : str
                base prefix, i.e. '/Desktop/SB_project1_prmatrix'
        """

        # construct savepath
        thedate = datetime.now().strftime("%y%m%d")
        savepath = filename + '_' + self.status + '_' + thedate + '.tsv'
        dated_filename = savepath.split('/')[-1]
        print(dated_filename)
        # save proteome as .tsv
        #self.proteome.to_csv(savepath, sep = '\t', index = True)
        return 
    
    def save_joblib():
        return
