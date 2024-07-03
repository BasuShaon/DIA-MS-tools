# Data-independent acquisition (DIA) tools

Automated pipeline for post-processing proteomes acquired with DIA-MS. Modules include QC, filtration, imputation, batch correction, and a wrapper for maxLFQ algorithm.

## Requirements

- Python 3.11 (anndata, datetime, gseapy, math, matplotlib, numpy, scanpy, scipy, seaborn, scikit-learn, statsmodels) 

- R 4.3.1 / Bioconductor 3.18 (AnnotationDbi, data.table, org.Hs.eg.db, iq)

## Installation 

Clone 'DIA-MS-tools' locally and cd into directory. Install requirements in 'requirements.txt' using pip / venv.

## Directory setup

Make 2 sub folders within cloned repo directory (named 'input' 'output') and copy following files into 'input':

### DIA-MS output (.tsv)

Precursor matrix containing sample IDs as columns and precursor intensities as rows. Rows should be indexed by protein and fragment ID (first and second columns). 

### Metadata (.tsv)

Metadata file which contains batch IDs as a categorical feature (column) and sample IDs as rows. Target column must be named 'MS.Batch'.

## Execution

Execute main.py from root of cloned directory.

python3 main.py -m <matrix_file> -b <batch_data_file> -o <output_prefix>

![output directory](https://github.com/BasuShaon/DIA-MS-tools/blob/main/screen.png)
