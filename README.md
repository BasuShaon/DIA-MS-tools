# DIA-MS-tools

![output directory](https://github.com/BasuShaon/DIA-MS-tools/blob/main/docs/pipe.png)

## Description

Data-driven pipeline to build complete proteomes with reduced batch effects. Resulting abundance follows a normal probability distribution. 

Designed for high-throughput proteomic drug screens in mammalian cells using data-independent-acquisition mass-spectrometry (DIA-MS). Tested thus far on fragments detected from 1000, 2000, and 5000+ sample drug screens (subjected to DAI-MS [dia-PASEF / scanning SWATH] + DIA-NN).

It should work using detected fragments from any DIA-MS setup in theory. More samples (1000+) & fragments (500000+) the better. 

Modules include QC, filtration, imputation, batch correction, and a wrapper for maxLFQ algorithm.

## Requirements

- Python 3.11.5 (download dependencies in local or virtual environment from 'requirements.txt') 

- R 4.3.1 (downloads required packages from CRAN at run time)

## Installation 

Clone 'DIA-MS-tools' locally and cd into root directory. Before running, make an input subfolder (named 'input'), and copy the following two files:

### MS data (.tsv)

Precursor matrix containing sample IDs as columns and precursor intensities as rows. Rows should be indexed by protein and fragment ID (first and second columns). 

### Batch data (.tsv)

Metadata file which contains batch IDs as a categorical feature (column) and sample IDs as rows. Target column must be named 'MS.Batch'.

## Execution

Create an output subfolder in directory (named 'output'), and execute main.py.

python3 main.py -m <matrix_file> -b <batch_data_file> -o <output_prefix>

![output directory](https://github.com/BasuShaon/DIA-MS-tools/blob/main/docs/screen.png)

