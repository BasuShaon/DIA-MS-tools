# %% Most imports for all the devices
from datetime import datetime
import pandas as pd
import numpy as np

# Device with global functions (CV calc, plotting, saving, ...) to inherit
class Carrier: 
    """ 
    Carrier class that contains information for processing proteomes 

    Attributes
    ----------
    frame : pd.DataFrame
        precursor intensity data 

    status: str
        status of the carrier frame
    
    batch_vec : pd.Series
        batch data

    save_frame(filename): fuction
        saves frame as .tsv

    """

    def __init__(self, filename, intensity_data, batch_data):

        self.frame = intensity_data
        self.status = 'raw'
        self.batch_vec = batch_data

    def save_frame(self, filename):
        """  
        Saves carrier frame as a .tsv with status and date suffix

        Parameters:
        --------       
        fname : str
            base filename, i.e. '/Desktop/SB_project1_prmatrix'

        """

        # construct savepath
        thedate = datetime.now().strftime("%y%m%d")
        savepath = filename + '_' + self.status + '_' + thedate + '.tsv'
        dated_filename = savepath.split('/')[-1]
        print(dated_filename)
        # save as carrier frame as .tsv
        self.frame.to_csv(savepath, sep = '\t', index = True)
        return 
    
