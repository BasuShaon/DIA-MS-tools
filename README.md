# DIA-MS-tools

Automated pipeline for pre-processing proteomes acquired with data-independent acquisition mass spectrometry (DIA-MS). Modules include QC, filtration, imputation, batch correction, and a wrapper for maxLFQ algorithm.

## Requirements

- Python 3.11.5 (install dependencies in local or virtual environment from 'requirements.txt')

- R 4.3.1 (installs required packages from CRAN at run time)

## Installation 

Clone 'DIA-MS-tools' locally and cd into the directory. Before running, make an input subfolder within the directory (named 'input'), and copy the following two files:

### DIA-MS output (.tsv)

Precursor matrix containing sample IDs as columns and precursor intensities as rows. Rows should be indexed by protein and fragment ID (first and second columns). 

### Metadata (.tsv)

Metadata file which contains batch IDs as a categorical feature (column) and sample IDs as rows. Target column must be named 'MS.Batch'.

## Execute python script (main.py)

Create an output subfolder in directory (named 'output'), and execute main.py.

python3 main.py -m <matrix_file> -b <batch_data_file> -o <output_prefix>

![output directory](https://github.com/BasuShaon/DIA-MS-tools/blob/main/screen.png)
