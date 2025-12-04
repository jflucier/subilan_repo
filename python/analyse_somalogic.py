import os.path
import subprocess
import argparse
import sys
import os
import re
import pandas as pd
import numpy as np
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns

def format_adat(adat_f,samples_f,comp_name,outpath):
    print("Formatting adat in R")
    command = (
            "Rscript --vanilla "
            + os.path.join(os.path.dirname(__file__), '../R/somalogic_to_tsv.R ')
            + " -o " + outpath
            + " -a " + adat_f
            + " -s " + samples_f
            + " -c " + comp_name
    )
    print(command)
    return_code = subprocess.call(command, shell=True)

    if return_code == 0:
        print("somalogic_to_tsv executed successfully.")
    else:
        print("somalogic_to_tsv failed with return code", return_code)
        exit(1)

    formatted_adat = pd.read_csv(os.path.join(outpath, 'adat.raw.t.tsv'), sep="\t")
    return formatted_adat

def gen_stats(adat, targets, comp_name, lbl1, lbl2, outpath):

    cols=[comp_name]
    seqid_cols = [col for col in adat.columns.tolist() if re.search(r"seq\..*",col)]
    cols.extend(seqid_cols)
    filtered_adat = adat[cols]

    selected_rows = filtered_adat[(filtered_adat[comp_name] == lbl1) | (filtered_adat[comp_name] == lbl2)]

    # Group by 'group' column and calculate the mean for each group
    group_means = selected_rows.groupby(comp_name).mean()
    # Corrected section for swapping
    group_means_transposed = group_means.T
    group_means_transposed = group_means_transposed.rename(columns={lbl1: f"avg_{lbl1}", lbl2: f"avg_{lbl2}"})

    # Group by 'group' column and calculate the standard deviation for each group
    group_std = selected_rows.groupby(comp_name).std()
    group_std_transposed = group_std.T
    group_std_transposed = group_std_transposed.rename(columns={lbl1: f"stdev_{lbl1}", lbl2: f"stdev_{lbl2}"})

    # Compute log2 fold change of the group means
    log2_fold_change = np.log2(group_means.loc[lbl1] / group_means.loc[lbl2])
    stats = log2_fold_change.reset_index()
    stats.columns = ['feature', 'log2fc']
    # group_means = group_means.reset_index()
    # group_means_transposed = group_means.T
    # group_means_transposed.columns = [f"avg_{lbl1}", f"avg_{lbl2}"]
    stats = stats.merge(group_means_transposed, left_on='feature', right_index=True, how='left')
    # group_std = group_std.reset_index()
    # group_std_transposed = group_std.T
    # group_std_transposed.columns = [f"stdev_{lbl1}", f"stdev_{lbl2}"]
    stats = stats.merge(group_std_transposed, left_on='feature', right_index=True, how='left')


    # Calculate t-test p-values for each feature
    group1_data = selected_rows[selected_rows[comp_name] == lbl1].drop(columns=comp_name)
    group2_data = selected_rows[selected_rows[comp_name] == lbl2].drop(columns=comp_name)
    p_values = []
    for col in group1_data.columns:
        t_stat, p_val = ttest_ind(group1_data[col], group2_data[col], nan_policy='omit')
        p_values.append(p_val)
    stats['p_value'] = p_values

    # Adjust p-values using Benjamini-Hochberg procedure
    adjusted_p_values = multipletests(p_values, method='fdr_bh')[1]
    stats['adjusted_p_value'] = adjusted_p_values


    stats_ann = targets.merge(stats, how='left', left_on='AptName', right_on='feature')

    stats_ann.to_csv(os.path.join(outpath, f'{lbl1}vs{lbl2}.stats.tsv'), sep='\t', index=False)
    return stats_ann

def gen_volcano(stats_ann, comp_name, lbl1, lbl2,outpath):

    stats_ann['-log10(adj.p_value)'] = -np.log10(stats_ann['adjusted_p_value'])
    conditions = [
        (stats_ann['log2fc'] <= -1) & (stats_ann['-log10(adj.p_value)'] >= 1.3),
        (stats_ann['log2fc'] >= 1) & (stats_ann['-log10(adj.p_value)'] >= 1.3)
    ]
    choices = ['Down', 'Up']
    stats_ann['color'] = np.select(conditions, choices, default='NS')

    # Create volcano plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=stats_ann, x='log2fc', y='-log10(adj.p_value)', hue='color', palette=['grey', 'blue', 'red'])

    for i in range(stats_ann.shape[0]):
        if stats_ann.loc[i, 'color'] == 'Down':
            plt.text(
                stats_ann.loc[i, 'log2fc'],
                stats_ann.loc[i, '-log10(adj.p_value)'],
                stats_ann.loc[i, 'EntrezGeneSymbol'],
                fontsize=9, ha='right'
            )
        if merged_df.loc[i, 'color'] == 'Up':
            plt.text(
                stats_ann.loc[i, 'log2fc'],
                stats_ann.loc[i, '-log10(adj.p_value)'],
                stats_ann.loc[i, 'EntrezGeneSymbol'],
                fontsize=9, ha='left'
            )

    # Add labels and title
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10 Adjusted p-value')
    plt.title(f"{lbl1} vs {lbl2}")
    plt.axhline(y=1.3, color='black', linestyle='--', linewidth=0.7)
    plt.axvline(x=-1, color='black', linestyle='--', linewidth=0.7)
    plt.axvline(x=1, color='black', linestyle='--', linewidth=0.7)

    # plt.show()
    plt.savefig(os.path.join(outpath, f'{lbl1}vs{lbl2}.volcano.pdf'))
    plt.close()


def gen_pca(adat, comp_name, lbl1, lbl2, outpath):
    # Perform PCA
    pca = PCA(n_components=2)

    cols=[comp_name]
    seqid_cols = [col for col in adat.columns.tolist() if re.search(r"seq\..*",col)]
    cols.extend(seqid_cols)
    filtered_adat = adat[cols]

    selected_rows = filtered_adat[(filtered_adat[comp_name] == lbl1) | (filtered_adat[comp_name] == lbl2)]

    selected_rows_data = selected_rows.drop(columns=comp_name)
    pca_result = pca.fit_transform(selected_rows_data)

    pca_df = pd.DataFrame(pca_result, columns=['PC1', 'PC2'])
    pca_df[comp_name] = selected_rows[comp_name].values

    # Create PCA plot
    plt.figure(figsize=(10, 6))
    sns.scatterplot(data=pca_df, x='PC1', y='PC2', hue=comp_name, palette=['blue', 'red'])

    # Add labels and title
    plt.xlabel(f'PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)')
    plt.ylabel(f'PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)')
    plt.title(f"PCA Plot: {lbl1} vs {lbl2}")

    # plt.show()
    plt.savefig(os.path.join(outpath, f'{lbl1}vs{lbl2}.pca.pdf'))
    plt.close()

def gen_heatmap(adat, stats_ann, comp_name, lbl1, lbl2, outpath):

    cols=[comp_name]
    seqid_cols = [col for col in adat.columns.tolist() if re.search(r"seq\..*",col)]
    cols.extend(seqid_cols)
    filtered_adat = adat[cols]

    selected_rows = filtered_adat[(filtered_adat[comp_name] == lbl1) | (filtered_adat[comp_name] == lbl2)]
    df_no_seq = selected_rows.filter(regex='^(?!seq_).*$')
    transposed_df = df_no_seq.T
    transposed_df.to_csv(os.path.join(outpath, f'{lbl1}vs{lbl2}.allrows.tsv'), sep='\t', index=False)

    # Filter features based on conditions
    filtered_features = stats_ann[
        (stats_ann['log2fc'] <= -1) & (stats_ann['adjusted_p_value'] <= 0.05) |
        (stats_ann['log2fc'] >= 1) & (stats_ann['adjusted_p_value'] <= 0.05)
    ]['feature'].tolist()

    heatmap_data = selected_rows[selected_rows.columns.intersection(filtered_features)].copy()
    heatmap_data[comp_name] = selected_rows[comp_name].values

    # Normalize the heatmap data (z-score normalization)
    heatmap_data_norm = heatmap_data.drop(columns=[comp_name]).apply(lambda x: (x - x.mean()) / x.std(), axis=0)
    heatmap_data_norm[comp_name] = heatmap_data[comp_name].values

    # Map the features to EntrezGeneSymbol for row labels
    feature_to_symbol = stats_ann.set_index('feature')['EntrezGeneSymbol'].to_dict()
    heatmap_data_norm.rename(columns=feature_to_symbol, inplace=True)

    # Create a heatmap
    heatmap_data_norm.set_index(comp_name, inplace=True)
    plt.figure(figsize=(10, 8))
    clustermap = sns.clustermap(heatmap_data_norm.T, cmap='coolwarm', annot=False)
    clustermap.ax_heatmap.collections[0].colorbar.set_label('Z-Score')

    # Adjust layout to provide more space for the title and scale bar
    clustermap.ax_heatmap.figure.subplots_adjust(top=0.92)
    clustermap.fig.suptitle(f"Heatmap: {lbl1} vs {lbl2}", x=0.5, ha='center')

    # Adjust color bar width
    cbar = clustermap.ax_heatmap.collections[0].colorbar
    cbar.ax.tick_params(labelsize=8)
    cbar.ax.set_aspect(2)

    # Save clustermap
    plt.savefig(os.path.join(outpath, f'{lbl1}vs{lbl2}.heatmap.pdf'))
    plt.close()

if __name__ == '__main__':

    argParser = argparse.ArgumentParser()

    # mandatory
    argParser.add_argument(
        "-i",
        "--adat",
        help="adat file provided by somalogic",
        required=True
    )

    argParser.add_argument(
        "-s",
        "--samples",
        help="Samples annotation file",
        required=True
    )

    argParser.add_argument(
        "-o",
        "--outpath",
        help="output path",
        required=True
    )

    argParser.add_argument(
        "-g",
        "--groups",
        help="Groups column name",
        required=True
    )

    argParser.add_argument(
        "-g1",
        "--group1",
        help="group1 name. Will be numerator in log2 fold calculation: log2(g1/g2)",
        required=True
    )

    argParser.add_argument(
        "-g2",
        "--group2",
        help="group2 name. Will be denominator in log2 fold calculation: log2(g1/g2)",
        required=True
    )

    args = argParser.parse_args()

    adat_f = args.adat
    samples_f = args.samples
    outpath = args.outpath
    groups_col = args.groups
    group1 = args.group1
    group2 = args.group2

    # adat_f = "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take3/SS-2453319_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.20240528.adat"
    # samples_f = "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take3/plate_sample_annot.tsv"
    # outpath = "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take3"
    # groups_col = "CompE"
    # group1 = "EPost"
    # group2 = "EPre"

    if not os.path.exists(outpath):
        os.makedirs(outpath, exist_ok=True, recursive=True)

    print("Formatting somalogic Adat")
    formated_adat = format_adat(
        adat_f,
        samples_f,
        groups_col,
        outpath
    )

    print("Reading targets")
    targets = pd.read_csv(os.path.join(outpath, 'all_targets.tsv'), sep="\t")

    print("Generate stats")
    stats_ann = gen_stats(
        formated_adat,
        targets,
        groups_col,
        group1,
        group2,
        outpath
    )

    print("Generate PCA plot")
    gen_pca(
        formated_adat,
        groups_col,
        group1,
        group2,
        outpath
    )

    print("Generate heatmap plot")
    gen_heatmap(
        formated_adat,
        stats_ann,
        groups_col,
        group1,
        group2,
        outpath
    )

