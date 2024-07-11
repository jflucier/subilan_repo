import argparse
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import numpy as np
# import scipy as sp
import scipy.stats as stats

def analyse_tcga(report, comp, ref, o):
    all_data = pd.read_csv(report, sep='\t', index_col="sample")
    comparisons = pd.read_csv(comp, sep='\t', index_col="label")
    res = None
    for label, row in comparisons.iterrows():
        print(f"running {label}")
        # filter data based on category field
        pd_rep = all_data.loc[(all_data['category'] == row['c_category']) | (all_data['category'] == row['n_category'])]
        # pd_rep = pd.read_csv(report, sep='\t', index_col="sample")
        low, high = pd_rep[ref].quantile([0.25, 0.75])
        #nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
        nd = pd_rep.loc[(pd_rep[ref] <= low) | (pd_rep[ref] >= high)]
        nd = nd.drop(columns=["category"])
        x = nd.corr(method='spearman')
        tmp = x['SOCS1'].to_frame()
        tmp = tmp.rename(columns={"SOCS1": label})
        if res is None:
            res = tmp
        else:
            res = res.merge(tmp, left_index=True, right_index=True, how='outer')
        # res.join(x['SOCS1'].set_axis(df1.index))
    res.to_csv(f"{o}/correlation/correlation.all.spearman.tsv", sep="\t")

def analyse_tcga_high(report, comp, ref, o):
    all_data = pd.read_csv(report, sep='\t', index_col="sample")
    comparisons = pd.read_csv(comp, sep='\t', index_col="label")
    res = None
    for label, row in comparisons.iterrows():
        print(f"running {label}")
        # filter data based on category field
        pd_rep = all_data.loc[(all_data['category'] == row['c_category']) | (all_data['category'] == row['n_category'])]
        # pd_rep = pd.read_csv(report, sep='\t', index_col="sample")
        low, high = pd_rep[ref].quantile([0.25, 0.75])
        #nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
        nd = pd_rep.loc[(pd_rep[ref] >= high)]
        nd = nd.drop(columns=["category"])
        x = nd.corr(method='spearman')
        tmp = x['SOCS1'].to_frame()
        tmp = tmp.rename(columns={"SOCS1": label})
        if res is None:
            res = tmp
        else:
            res = res.merge(tmp, left_index=True, right_index=True, how='outer')
        # res.join(x['SOCS1'].set_axis(df1.index))
    res.to_csv(f"{o}/correlation/correlation.high.spearman.tsv", sep="\t")

def analyse_tcga_low(report, comp, ref, o):
    all_data = pd.read_csv(report, sep='\t', index_col="sample")
    comparisons = pd.read_csv(comp, sep='\t', index_col="label")
    res = None
    for label, row in comparisons.iterrows():
        print(f"running {label}")
        # filter data based on category field
        pd_rep = all_data.loc[(all_data['category'] == row['c_category']) | (all_data['category'] == row['n_category'])]
        # pd_rep = pd.read_csv(report, sep='\t', index_col="sample")
        low, high = pd_rep[ref].quantile([0.25, 0.75])
        #nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
        nd = pd_rep.loc[(pd_rep[ref] <= low)]
        nd = nd.drop(columns=["category"])
        x = nd.corr(method='spearman')
        tmp = x['SOCS1'].to_frame()
        tmp = tmp.rename(columns={"SOCS1": label})
        if res is None:
            res = tmp
        else:
            res = res.merge(tmp, left_index=True, right_index=True, how='outer')
    res.to_csv(f"{o}/correlation/correlation.low.spearman.tsv", sep="\t")

def make_boxplot(report, comp, ref, o):
    all_data = pd.read_csv(report, sep='\t', index_col="sample")
    socs1_data = all_data.filter(items=['category', ref])
    comparisons = pd.read_csv(comp, sep='\t', index_col="label")
    for label, row in comparisons.iterrows():
        print(f"running {label}")
        # filter data based on category field
        pd_rep_c = socs1_data.loc[
            (socs1_data['category'] == row['c_category'])]
        pd_rep_n = socs1_data.loc[
            (socs1_data['category'] == row['n_category'])]
        pval = stats.mannwhitneyu(x=pd_rep_c[ref], y=pd_rep_n[ref], alternative='two-sided')
        comparisons.loc[label, 'MW_pval'] = pval.pvalue
        comparisons.loc[label, 'N_mean'] = pd_rep_n[ref].mean()
        comparisons.loc[label, 'C_mean'] = pd_rep_c[ref].mean()
        comparisons.loc[label, 'C-N'] = comparisons.loc[label, 'C_mean'] - comparisons.loc[label, 'N_mean']

        pd_rep = socs1_data.loc[(socs1_data['category'] == row['c_category']) | (socs1_data['category'] == row['n_category'])]
        myFig = plt.figure()
        ax = myFig.add_subplot(111)
        ax = pd_rep.boxplot(by='category', ax=ax)
        ax.get_figure().suptitle(f"{label} Cancer (N={row['c_samplenbr']}) vs Normal (N={row['n_samplenbr']}). pval={pval.pvalue:.2e}")
        myFig.savefig(f"{o}/boxplot/tcga_{ref}_{label}.png",
                      format="png", bbox_inches="tight")
        plt.close()

    comparisons.to_csv(f"{o}/boxplot/tcga_{ref}_stats.tsv", sep="\t")

def analyse_tcga_cancer_only(report, comp, ref, o):
    all_data = pd.read_csv(report, sep='\t', index_col="sample")
    comparisons = pd.read_csv(comp, sep='\t', index_col="label")
    for label, row in comparisons.iterrows():
        print(f"running {label}")
        # filter data based on category field
        pd_rep = all_data.loc[(all_data['category'] == row['c_category'])]
        # pd_rep = pd.read_csv(report, sep='\t', index_col="sample")
        low, high = pd_rep[ref].quantile([0.25, 0.75])
        # nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
        nd = pd_rep.loc[(pd_rep[ref] <= low) | (pd_rep[ref] >= high)]
        nd = nd.drop(columns=["category"])
        x = nd.corr(method='spearman')
        x.to_csv(f"{o}/{label}.cancer.spearman.tsv", sep="\t")


def analyse_tcga_normal_only(report, comp, ref, o):
    all_data = pd.read_csv(report, sep='\t', index_col="sample")
    comparisons = pd.read_csv(comp, sep='\t', index_col="label")
    for label, row in comparisons.iterrows():
        print(f"running {label}")
        # filter data based on category field
        pd_rep = all_data.loc[(all_data['category'] == row['n_category'])]
        # pd_rep = pd.read_csv(report, sep='\t', index_col="sample")
        low, high = pd_rep[ref].quantile([0.25, 0.75])
        # nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
        nd = pd_rep.loc[(pd_rep[ref] <= low) | (pd_rep[ref] >= high)]
        nd = nd.drop(columns=["category"])
        x = nd.corr(method='spearman')
        x.to_csv(f"{o}/{label}.normal.spearman.tsv", sep="\t")


if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument(
        "-ug",
        "--diann_ug_report",
        help="DIANN unqiue genen report",
        required=True
    )

    argParser.add_argument(
        "-s",
        "--samplesheet",
        help="TSV file with these columns: name<tab>group<tab>file",
        required=True
    )

    argParser.add_argument(
        "-o",
        "--out_prefix",
        help="output path prefix",
        required=True
    )

    args = argParser.parse_args()

    analyse_diann_ug(
        args.diann_ug_report,
        args.samplesheet,
        args.reference,
        args.out_prefix
    )

    analyse_tcga_low(
        args.tcga_report,
        args.comparisons,
        args.reference,
        args.out_prefix
    )

    analyse_tcga(
        args.tcga_report,
        args.comparisons,
        args.reference,
        args.out_prefix
    )