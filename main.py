"""
Proteomics Data Processing Pipeline
"""
import sys
import os
import pandas as pd
import argparse

# Add src to path
sys.path.append(os.path.abspath("code"))

from carrier import Carrier
import detectioncontrol
import imputer
import combat
import precursor2protein

def main(data_df, meta_df, proj, alpha, bound):
    # Configuration
    DATA_DF = data_df
    META_DF = meta_df
    PROJECT_NAME = proj
    ALPHA = alpha
    BOUND = bound
    DATA_DIR = os.path.join(os.path.dirname(__file__), 'input')
    OUT_DIR = os.path.join(os.path.dirname(__file__), 'output')

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
    
    print("3. Running detection control...")
    optimum = detectioncontrol.calculate_optimum_threshold(shogoki, alpha=ALPHA)
    detectioncontrol.detection_control(shogoki, optimum)
    shogoki.save()
    
    print("4. Running imputation...")
    imputer.preprocess_data(shogoki)
    imputer.compute_log2_means_and_missingness(shogoki)
    if bound == 'Auto':
        bound, glm = imputer.detection_probability_curve(shogoki)
        imputer.mixed_imputation_in_batch(shogoki, bound)
    else:
        bound == BOUND
    shogoki.save()
    
    print("5. Running batch correction...")
    combat.batch_correct(shogoki)
    combat.CV_plots(shogoki, shogoki.proteome_log2_beforecombat, 'before Combat')
    combat.CV_plots(shogoki, shogoki.proteome, 'after Combat')
    shogoki.save()
    
    # # MaxLFQ
    # print("6. Running MaxLFQ...")
    # precursor2protein.format_df_maxlfq(shogoki)
    # shogoki.save()
    
    # output_path = os.path.join(
    #     shogoki.outer_path,
    #     f"{shogoki.projectname}_{shogoki.status}_{shogoki.thedate}.tsv"
    # )
    # precursor2protein.maxlfq(shogoki, output_path)
    
    # print(f"\n✓ Pipeline complete! Status: {shogoki.status}")


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
    parser.add_argument(
        "-a", type=float, default=0.99,
        help="Alpha/significance parameter for detection control (default: 0.99)."
    )
    parser.add_argument(
        "-b", default="Auto",
        help="Boundary for MNAR/MAR split. Float or 'Auto' to fit via logistic curve (default: Auto)."
    )

    args = parser.parse_args()
    main(
        data_df=args.d,
        meta_df=args.m,
        proj=args.p,
        alpha=args.a,
        bound=args.b,
    )