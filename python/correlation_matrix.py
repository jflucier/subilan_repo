import argparse
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import numpy as np
# import scipy as sp
import scipy.stats as stats
import seaborn as sns

def gen_corr_matrix(report, o):
    all_data = pd.read_csv(report, sep='\t', index_col="Genes")

    all_data_t = all_data.transpose()
    corrM = all_data.corr()

    x = sns.heatmap(corrM, cmap="Blues", annot=True )
    x.figure.savefig(o)

if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument(
        "-ug",
        "--diann_ug_report",
        help="DIANN unqiue gene report",
        required=True
    )

    argParser.add_argument(
        "-o",
        "--out_prefix",
        help="output path prefix",
        required=True
    )

    args = argParser.parse_args()

    gen_corr_matrix(
        args.diann_ug_report,
        args.out_prefix
    )