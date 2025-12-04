import argparse
import matplotlib.pyplot as plt
import pandas as pd
from pandas import DataFrame
import numpy as np
# import scipy as sp
import scipy.stats as stats

def estimates_filter(pd_rep,label):
    # A. Create a temporary 'group' column based on the index (sample ID: 'TCGA-...' or 'GTEX-...')
    # We assume the group name is the part before the first hyphen in the sample ID.
    pd_rep['group'] = pd_rep.index.str.split('-').str[0]

    cutoffs = pd_rep.groupby('group')['ESTIMATE_score'].quantile(0.75)

    # C. Print the calculated cutoffs to the console
    print(f"\nCutoffs for comparison '{label}':")
    print(cutoffs)

    # B. Group the data by 'group' and apply the Q4 filter (ESTIMATE_score <= Q3)
    # group_keys=False prevents the temporary 'group' column from becoming part of the index.
    pd_rep = pd_rep.groupby('group', group_keys=False).apply(
        lambda x: x[x['ESTIMATE_score'] <= x['ESTIMATE_score'].quantile(0.75)],
        include_groups=False  # <-- ADD THIS LINE
    )
    return pd_rep

def run_correlation(label, pd_rep, ref):
    low, high = pd_rep[ref].quantile([0.25, 0.75])
    # nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
    nd = pd_rep.loc[(pd_rep[ref] <= low) | (pd_rep[ref] >= high)]
    nd = nd.drop(columns=["category"])
    x = nd.corr(method='spearman')
    tmp = x[ref].to_frame()
    tmp = tmp.rename(columns={ref: label})
    return tmp

def run_splicing(o, label, pd_rep, splicing_data, ref):
    low, high = pd_rep[ref].quantile([0.25, 0.75])
    # Filter 'pd_rep' to get sample IDs for low- and high-expression groups
    low_samples = pd_rep.loc[pd_rep[ref] <= low].index.tolist()
    high_samples = pd_rep.loc[pd_rep[ref] >= high].index.tolist()

    # Combine sample lists for filtering the large splicing data
    all_filtered_samples = low_samples + high_samples
    # Ensure all filtered samples exist in the splicing data columns to avoid errors
    existing_samples = [s for s in all_filtered_samples if s in splicing_data.columns]
    # The result 'pd_splice' has Transcript IDs as index and Sample IDs as columns.
    pd_splice = splicing_data[existing_samples].copy()

    # Samples belonging to the high expression group
    high_splice = pd_splice[[s for s in high_samples if s in pd_splice.columns]]
    # Samples belonging to the low expression group
    low_splice = pd_splice[[s for s in low_samples if s in pd_splice.columns]]

    # 1. Transpose the dataframes so samples are rows (as required for output)
    high_t = high_splice.T
    low_t = low_splice.T

    # 2. Add the reference group label column to each set
    high_t['ref_group'] = 'HIGH'
    low_t['ref_group'] = 'LOW'

    # 3. Concatenate the two dataframes into a single result
    combined_raw_data = pd.concat([high_t, low_t])
    combined_raw_data = combined_raw_data.sort_values(by='ref_group', ascending=False)
    combined_raw_data.to_csv(f"{o}/tcga.splicing.{label}.raw.tsv", sep='\t', index=True)

    # Calculate the row-wise (transcript-wise) average expression
    high_avg = high_splice.mean(axis=1)
    low_avg = low_splice.mean(axis=1)

    transcript_difference = (high_avg - low_avg).to_frame(name='avg_diff')
    # Add a column for the label (e.g., "KIRC_Normal_vs_Tumor")
    transcript_difference[label] = transcript_difference['avg_diff']
    # Get the absolute values of the transcript differences
    abs_diff = transcript_difference['avg_diff'].abs()

    # Calculate the mean and standard deviation (stdev) of the absolute differences
    mean_abs_diff = abs_diff.mean()
    std_abs_diff = abs_diff.std()

    # Calculate the filtering threshold: mean + 1 * stdev
    # This filters out transcripts whose absolute change is below the majority.
    threshold = mean_abs_diff + (3 * std_abs_diff)

    print(f"  -> Splicing Delta-PSI Cutoff ({label}): mean(|diff|) + 1*stdev = {threshold:.4f}")
    filtered_difference = transcript_difference.loc[abs_diff >= threshold]
    filtered_difference = filtered_difference.sort_values(by='avg_diff', ascending=True)
    filtered_difference.to_csv(f"{o}/tcga.splicing.{label}.tsv", sep="\t")


def analyse_tcga(report, splicing_rep, comp, ref, o):
    all_data = pd.read_csv(report, sep='\t', index_col="sample")
    splicing_data = pd.read_csv(splicing_rep, sep='\t', index_col="sample")
    comparisons = pd.read_csv(comp, sep='\t', index_col="label")
    res_all = None
    res_tcga = None
    res_gtex = None
    for label, row in comparisons.iterrows():
        print(f"running {label} on all samples")
        # filter data based on category field using normal and cancer samples
        # pd_rep = all_data.loc[(all_data['category'] == row['c_category']) | (all_data['category'] == row['n_category'])].copy()
        # pd_rep = estimates_filter(pd_rep,label)
        # tmp = run_correlation(label, pd_rep, ref)
        # if res_all is None:
        #     res_all = tmp
        # else:
        #     res_all = res_all.merge(tmp, left_index=True, right_index=True, how='outer')

        # filter data based on category field using cancer samples
        print(f"running {label} on tcga samples")
        pd_rep = all_data.loc[(all_data['category'] == row['c_category'])].copy()
        pd_rep = estimates_filter(pd_rep,label)
        tmp = run_correlation(label, pd_rep, ref)
        if res_tcga is None:
            res_tcga = tmp
        else:
            res_tcga = res_tcga.merge(tmp, left_index=True, right_index=True, how='outer')

        run_splicing(o, label, pd_rep, splicing_data, ref)

        # filter data based on category field using gtex samples
        # print(f"running {label} on gtex samples")
        # pd_rep = all_data.loc[(all_data['category'] == row['n_category'])].copy()
        # pd_rep = estimates_filter(pd_rep,label)
        # print(f"\nRunning correlation for '{label}':")
        # tmp = run_correlation(label, pd_rep, ref)
        # if res_gtex is None:
        #     res_gtex = tmp
        # else:
        #     res_gtex = res_gtex.merge(tmp, left_index=True, right_index=True, how='outer')

    # res_all.to_csv(f"{o}/all.correlation.spearman.estimateQ4.tsv", sep="\t")
    res_tcga.to_csv(f"{o}/tcga.correlation.spearman.estimateQ4.tsv", sep="\t")
    # res_gtex.to_csv(f"{o}/gtex.correlation.spearman.estimateQ4.tsv", sep="\t")

# def analyse_tcga_high(report, comp, ref, o):
#     all_data = pd.read_csv(report, sep='\t', index_col="sample")
#     comparisons = pd.read_csv(comp, sep='\t', index_col="label")
#     res = None
#     for label, row in comparisons.iterrows():
#         print(f"running {label}")
#         # filter data based on category field
#         pd_rep = all_data.loc[(all_data['category'] == row['c_category']) | (all_data['category'] == row['n_category'])]
#         # pd_rep = pd.read_csv(report, sep='\t', index_col="sample")
#         low, high = pd_rep[ref].quantile([0.25, 0.75])
#         #nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
#         nd = pd_rep.loc[(pd_rep[ref] >= high)]
#         nd = nd.drop(columns=["category"])
#         x = nd.corr(method='spearman')
#         tmp = x['SOCS1'].to_frame()
#         tmp = tmp.rename(columns={"SOCS1": label})
#         if res is None:
#             res = tmp
#         else:
#             res = res.merge(tmp, left_index=True, right_index=True, how='outer')
#         # res.join(x['SOCS1'].set_axis(df1.index))
#     res.to_csv(f"{o}/correlation/correlation.high.spearman.tsv", sep="\t")
#
# def analyse_tcga_low(report, comp, ref, o):
#     all_data = pd.read_csv(report, sep='\t', index_col="sample")
#     comparisons = pd.read_csv(comp, sep='\t', index_col="label")
#     res = None
#     for label, row in comparisons.iterrows():
#         print(f"running {label}")
#         # filter data based on category field
#         pd_rep = all_data.loc[(all_data['category'] == row['c_category']) | (all_data['category'] == row['n_category'])]
#         # pd_rep = pd.read_csv(report, sep='\t', index_col="sample")
#         low, high = pd_rep[ref].quantile([0.25, 0.75])
#         #nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
#         nd = pd_rep.loc[(pd_rep[ref] <= low)]
#         nd = nd.drop(columns=["category"])
#         x = nd.corr(method='spearman')
#         tmp = x['SOCS1'].to_frame()
#         tmp = tmp.rename(columns={"SOCS1": label})
#         if res is None:
#             res = tmp
#         else:
#             res = res.merge(tmp, left_index=True, right_index=True, how='outer')
#     res.to_csv(f"{o}/correlation/correlation.low.spearman.tsv", sep="\t")
#
# def make_boxplot(report, comp, ref, o):
#     all_data = pd.read_csv(report, sep='\t', index_col="sample")
#     socs1_data = all_data.filter(items=['category', ref])
#     comparisons = pd.read_csv(comp, sep='\t', index_col="label")
#     for label, row in comparisons.iterrows():
#         print(f"running {label}")
#         # filter data based on category field
#         pd_rep_c = socs1_data.loc[
#             (socs1_data['category'] == row['c_category'])]
#         pd_rep_n = socs1_data.loc[
#             (socs1_data['category'] == row['n_category'])]
#         pval = stats.mannwhitneyu(x=pd_rep_c[ref], y=pd_rep_n[ref], alternative='two-sided')
#         comparisons.loc[label, 'MW_pval'] = pval.pvalue
#         comparisons.loc[label, 'N_mean'] = pd_rep_n[ref].mean()
#         comparisons.loc[label, 'C_mean'] = pd_rep_c[ref].mean()
#         comparisons.loc[label, 'C-N'] = comparisons.loc[label, 'C_mean'] - comparisons.loc[label, 'N_mean']
#
#         pd_rep = socs1_data.loc[(socs1_data['category'] == row['c_category']) | (socs1_data['category'] == row['n_category'])]
#         myFig = plt.figure()
#         ax = myFig.add_subplot(111)
#         ax = pd_rep.boxplot(by='category', ax=ax)
#         ax.get_figure().suptitle(f"{label} Cancer (N={row['c_samplenbr']}) vs Normal (N={row['n_samplenbr']}). pval={pval.pvalue:.2e}")
#         myFig.savefig(f"{o}/boxplot/tcga_{ref}_{label}.png",
#                       format="png", bbox_inches="tight")
#         plt.close()
#
#     comparisons.to_csv(f"{o}/boxplot/tcga_{ref}_stats.tsv", sep="\t")
#
# def analyse_tcga_cancer_only(report, comp, ref, o):
#     all_data = pd.read_csv(report, sep='\t', index_col="sample")
#     comparisons = pd.read_csv(comp, sep='\t', index_col="label")
#     for label, row in comparisons.iterrows():
#         print(f"running {label}")
#         # filter data based on category field
#         pd_rep = all_data.loc[(all_data['category'] == row['c_category'])]
#         # pd_rep = pd.read_csv(report, sep='\t', index_col="sample")
#         low, high = pd_rep[ref].quantile([0.25, 0.75])
#         # nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
#         nd = pd_rep.loc[(pd_rep[ref] <= low) | (pd_rep[ref] >= high)]
#         nd = nd.drop(columns=["category"])
#         x = nd.corr(method='spearman')
#         x.to_csv(f"{o}/{label}.cancer.spearman.tsv", sep="\t")
#
#
# def analyse_tcga_normal_only(report, comp, ref, o):
#     all_data = pd.read_csv(report, sep='\t', index_col="sample")
#     comparisons = pd.read_csv(comp, sep='\t', index_col="label")
#     for label, row in comparisons.iterrows():
#         print(f"running {label}")
#         # filter data based on category field
#         pd_rep = all_data.loc[(all_data['category'] == row['n_category'])]
#         # pd_rep = pd.read_csv(report, sep='\t', index_col="sample")
#         low, high = pd_rep[ref].quantile([0.25, 0.75])
#         # nd = pd_rep.query(f'{low}>{ref}' | f'{high}<{ref}')
#         nd = pd_rep.loc[(pd_rep[ref] <= low) | (pd_rep[ref] >= high)]
#         nd = nd.drop(columns=["category"])
#         x = nd.corr(method='spearman')
#         x.to_csv(f"{o}/{label}.normal.spearman.tsv", sep="\t")
#

if __name__ == '__main__':
    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument(
        "-r",
        "--tcga_report",
        help="tcga report filtered by category",
        required=True
    )

    argParser.add_argument(
        "-s",
        "--gtex_splicing_report",
        help="gtex splcing file TcgaTargetGtex_rsem_isopct",
        required=True
    )

    argParser.add_argument(
        "-comp",
        "--comparisons",
        help="category comparisons to perform",
        required=True
    )

    argParser.add_argument(
        "-ref",
        "--reference",
        help="reference column name",
        required=True
    )

    argParser.add_argument(
        "-o",
        "--out_prefix",
        help="output path prefix",
        required=True
    )

    args = argParser.parse_args()


    #
    # analyse_tcga_cancer_only(
    #     args.tcga_report,
    #     args.comparisons,
    #     args.reference,
    #     args.out_prefix
    # )
    #
    # analyse_tcga_normal_only(
    #     args.tcga_report,
    #     args.comparisons,
    #     args.reference,
    #     args.out_prefix
    # )

    # make_boxplot(
    #     args.tcga_report,
    #     args.comparisons,
    #     args.reference,
    #     args.out_prefix
    # )

    # analyse_tcga_high(
    #     args.tcga_report,
    #     args.comparisons,
    #     args.reference,
    #     args.out_prefix
    # )
    #
    # analyse_tcga_low(
    #     args.tcga_report,
    #     args.comparisons,
    #     args.reference,
    #     args.out_prefix
    # )

    analyse_tcga(
        args.tcga_report,
        args.gtex_splicing_report,
        args.comparisons,
        args.reference,
        args.out_prefix
    )