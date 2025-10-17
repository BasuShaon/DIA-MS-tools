# %% inherit device from abs path
import sys
import os
sys.path.append(os.path.join(os.path.dirname(__file__)))
from carrier import Carrier
import pickle
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

    wd = carrier.outerpath
    print(wd)
    
    file_path = longform
    carrier.status = carrier.status.removesuffix('_long') +  '_summarized'
    lfq_savename = carrier.projectname + '_' + carrier.status + '_' + carrier.thedate + '.tsv'
    print(lfq_savename)

    # Ensure convert is passed as a string "TRUE" or "FALSE"
    convert_str = "TRUE" if convert else "FALSE"
    command = ['Rscript', r_script_path, wd, file_path, lfq_savename, convert_str]
    process = subprocess.run(command, capture_output=True, text=True)
    print('STDOUT:', process.stdout)
    print('STDERR:', process.stderr)

    return


def main():
    # Configuration
    
    with open(os.path.join(os.path.dirname(__file__), "../output", "yongoki.pkl"), "rb") as f:
        gogoki = pickle.load(f)

    # MaxLFQ
    print("6. Running MaxLFQ...")
    format_df_maxlfq(gogoki)
    gogoki.save()
    
    output_path = os.path.join(
        gogoki.outerpath,
        f"{gogoki.projectname}_{gogoki.status}_{gogoki.thedate}.tsv"
    )
    
    maxlfq(gogoki, output_path)

    with open(os.path.join(
        gogoki.outerpath, 'gogoki.pkl'
    ),'wb') as f: 
        pickle.dump(gogoki, f)
    
if __name__ == "__main__":
    main()