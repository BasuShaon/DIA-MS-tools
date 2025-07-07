# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
from carrier import Carrier

# third party imports
import subprocess


# %% 
def format_df_maxlfq(carrier):

        corrected_proteome = carrier.proteome.T.reset_index().melt(id_vars=['Protein.Group', 'Precursor.Id'], 
                                var_name='Sample', 
                                value_name='Intensity')

        long_proteome = corrected_proteome.copy()

        long_proteome.columns = ['protein_list', 'id', 'sample_list', 'quant']

        long_proteome.set_index('id', inplace=True)

        long_proteome = long_proteome[['sample_list', 'protein_list', 'quant']]

        carrier.proteome = long_proteome

        carrier.status = carrier.status + '_long'

        return carrier

def maxlfq(carrier, longform, convert = False):

    r_script_path = os.path.join(os.path.dirname(__file__), 'maxLFQ.R')

    print(r_script_path)

    wd = carrier.outer_path

    print(wd)

    file_path = longform

    carrier.status = carrier.status + '_summarized'

    # Ensure convert is passed as a string "TRUE" or "FALSE"
    convert_str = "TRUE" if convert else "FALSE"

    command = ['Rscript', r_script_path, wd, file_path, carrier.status, convert_str]

    process = subprocess.run(command, capture_output=True, text=True)

    print('STDOUT:', process.stdout)

    print('STDERR:', process.stderr)

    return
