### https://www.sc-best-practices.org/preamble.html
### https://scanpy.readthedocs.io/en/stable/index.html
### https://github.com/scverse/scanpy

# import warnings
# warnings.filterwarnings("ignore", category=DeprecationWarning)

import os.path
import subprocess
import sys
import csv
import numpy as np
import scanpy as sc
import anndata as ad
import pandas as pd
import seaborn as sns

from matplotlib import pyplot as plt
from matplotlib.pyplot import rc_context
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import median_abs_deviation
from scipy import sparse, io
from scipy.sparse import csr_matrix, issparse

# import numba
# from numba.core.errors import NumbaDeprecationWarning, NumbaPendingDeprecationWarning
# warnings.simplefilter("ignore", category=NumbaDeprecationWarning)
import celltypist
from celltypist import models


sc.settings.set_figure_params(dpi=80, facecolor="white")
outpath = "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out"

def read_cellranger(samples):
    print("## Reading cellranger count matrix")
    adatas = {}
    for sample_id, filename in samples.items():
        print(f"Reading cellranger count matrix {filename}")
        sample_adata = sc.read_10x_h5(filename)
        print("Make unique gene names")
        sample_adata.var_names_make_unique()
        adatas[sample_id] = sample_adata
        print("Done")

    adata = ad.concat(adatas, label="sample")
    print("Make unique observation names")
    adata.obs_names_make_unique()
    print("Done")
    return adata

def output_qc(adata, output_pdf):

    if not os.path.exists(outpath):
        os.makedirs(outpath, exist_ok=True, recursive=True)

    with PdfPages(output_pdf) as pdf:
        # Inspect violin plots of some of the computed QC metrics:
        # 1 - the number of genes expressed in the count matrix
        # 2 - the total counts per cell
        # 3 - the percentage of counts in mitochondrial genes
        sc.pl.violin(
            adata,
            ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
            jitter=0.4,
            multi_panel=True,
            show=False
        )
        pdf.savefig(bbox_inches="tight")

        # scatter plot colored by pct_counts_mt
        sc.pl.scatter(
            adata,
            "total_counts",
            "n_genes_by_counts",
            color="pct_counts_mt",
            show=False
        )
        pdf.savefig(bbox_inches="tight")

        sc.pl.scatter(
            adata,
            "total_counts",
            "n_genes_by_counts",
            color="pct_counts_ribo",
            show=False
        )
        pdf.savefig(bbox_inches="tight")
        plt.close()

def output_feature_selection_qc(adata, fs_path):

    print("Outputting feature selection plots")

    with PdfPages(os.path.join(fs_path, "FeatureSelectionQC.pdf")) as pdf:

        fig, axes = plt.subplots(1, 3, figsize=(10, 5))
        p1 = sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axes[0])
        axes[0].set_title("Total counts")
        p2 = sns.histplot(adata.layers["log1p_norm"].sum(1), bins=100, kde=False, ax=axes[1])
        axes[1].set_title("Shifted logarithm")
        p3 = sns.histplot(
            adata.layers["analytic_pearson_residuals"].sum(1), bins=100, kde=False, ax=axes[2]
        )
        axes[2].set_title("Analytic Pearson residuals")
        pdf.savefig(bbox_inches="tight")

        # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # sc.pp.highly_variable_genes(adata, layer="analytic_pearson_residuals", n_top_genes=1000)
        # ax = sns.scatterplot(
        #     data=adata.var, x="means", y="dispersions", hue="highly_variable", s=5
        # )
        # # ax.set_xlim(None, 1.5)
        # # ax.set_ylim(-4, 4)
        # ax.set_xscale("log")
        # ax.set_yscale("log")
        # ax.set_title("Feature selection using Pearson residuals normalisation (from sc.experimental.pp.normalize_pearson_residuals)")
        # pdf.savefig(bbox_inches="tight")
        #

        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        ax = sns.scatterplot(
            data=adata.var, x="means", y="dispersions", hue="highly_variable", s=5
        )
        # ax.set_xlim(None, 1.5)
        # ax.set_ylim(-4, 4)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title("Feature selection using Pearson residuals normalisation (from highly_variable_genes flavor pearson_residuals)")
        pdf.savefig(bbox_inches="tight")

        fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        ax = sns.scatterplot(
            data=adata.var, x="means", y="dispersions", hue="highly_deviant", s=5
        )
        # ax.set_xlim(None, 1.5)
        # ax.set_ylim(-4, 4)
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_title(
            "Feature selection using Pearson residuals normalisation (from scry)")
        pdf.savefig(bbox_inches="tight")

        # fig, ax = plt.subplots(1, 1, figsize=(10, 5))
        # ax = sns.scatterplot(
        #     data=adata.var, x="means", y="dispersions", hue="highly_variable", s=5
        # )
        # # ax.set_xlim(None, 1.5)
        # # ax.set_ylim(-4, 4)
        # ax.set_xscale("log")
        # ax.set_yscale("log")
        # ax.set_title("Analytic Shifted logarithm normalisation")
        # pdf.savefig(bbox_inches="tight")

        plt.close()


def is_outlier(adata, metric: str, nmads: int):
    M = adata.obs[metric]
    outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
        np.median(M) + nmads * median_abs_deviation(M) < M
    )
    return outlier

def run_soupx(adata, soupx_path):
    # soupx_path = os.path.join(outpath, 'SoupX')
    # if not os.path.exists(soupx_path):
    #     os.makedirs(soupx_path, exist_ok=True)

    print(f"Preparing SoupX input files")
    adata_pp = adata.copy()
    sc.pp.normalize_per_cell(adata_pp)
    sc.pp.log1p(adata_pp)
    sc.pp.pca(adata_pp)
    sc.pp.neighbors(adata_pp)
    sc.tl.leiden(adata_pp, key_added="soupx_groups", flavor="igraph", n_iterations=2)

    # Preprocess variables for SoupX
    soupx_groups = adata_pp.obs["soupx_groups"]
    del adata_pp
    cells = adata.obs_names
    genes = adata.var_names
    data = adata.X.T

    print("Reading cellranger raw count matrices")
    adatas_raw = {}
    for sample_id, filename in samples_raw.items():
        print(f"Reading cellranger raw count matrix {filename}")
        adata_raw = sc.read_10x_h5(filename)
        print("Make unique gene names")
        adata_raw.var_names_make_unique()

        # del adata_raw
        adatas_raw[sample_id] = adata_raw
        print("Done")

    adatas_raw = ad.concat(adatas_raw, label="sample")
    print("Make raw data unique observation names")
    adatas_raw.obs_names_make_unique()
    data_tod = adatas_raw.X.T

    print("Saving filtered data matrix")
    sparse.save_npz(os.path.join(soupx_path, 'data_matrix.npz'), data)
    pd.Series(adata.obs_names.to_series()).to_csv(
        os.path.join(soupx_path, 'data_barcodes.tsv'),
        sep="\t", index=False, header=False
    )
    pd.Series(adata.var_names.to_series()).to_csv(
        os.path.join(soupx_path, 'data_features.tsv'),
        sep="\t", index=False,header=False
    )
    soupx_groups.to_csv(os.path.join(soupx_path, 'clusters.tsv'), sep="\t", index=True, header=False)

    print("Saving raw data matrix")
    sparse.save_npz(os.path.join(soupx_path, 'data_tod_matrix.npz'), data_tod)

    print("Running SoupX")
    command = (
            "Rscript --vanilla "
            + os.path.join(os.path.dirname(__file__), '../R/run_soupx.R ')
            + " -o " + os.path.join(soupx_path, 'soupx_matrix.npz')
            + " -f " + os.path.join(soupx_path, 'data_matrix.npz')
            + " -r " + os.path.join(soupx_path, 'data_tod_matrix.npz')
            + " -c " + os.path.join(soupx_path, 'clusters.tsv')
            + " -g " + os.path.join(soupx_path, 'data_features.tsv')
            + " -b " + os.path.join(soupx_path, 'data_barcodes.tsv')
            + " -x " + sys.executable
    )
    print(command)
    return_code = subprocess.call(command, shell=True)

    if return_code == 0:
        print("SoupX executed successfully.")
    else:
        print("SoupX failed with return code", return_code)

    print("Loading SoupX result matrix")
    soupx_data = sparse.load_npz(os.path.join(soupx_path, 'soupx_matrix.npz'))
    adata.layers["counts"] = adata.X
    adata.layers["soupX_counts"] = soupx_data.transpose()
    adata.X = adata.layers["soupX_counts"]
    return adata

def run_qc(adata, outpath):
    ## Quality Control
    print("## Running QC")
    qc_path = os.path.join(outpath, '01-QC')
    if not os.path.exists(qc_path):
        os.makedirs(qc_path, exist_ok=True)

    # mitochondrial genes, "MT-" for human, "Mt-" for mouse
    adata.var["mt"] = adata.var_names.str.startswith("mt-")
    # ribosomal genes
    adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True, percent_top=[20]
    )

    print("QC graph before filtering")
    output_qc(
        adata,
        os.path.join(qc_path, '01-BeforeFiltering.pdf')
    )

    print("Outlier filtering using 5 x MAD for log1p_total_counts, log1p_n_genes_by_counts and pct_counts_in_top_20_genes")

    adata.obs["outlier"] = (
            is_outlier(adata, "log1p_total_counts", 5)
            | is_outlier(adata, "log1p_n_genes_by_counts", 5)
            | is_outlier(adata, "pct_counts_in_top_20_genes", 5)
    )
    print(f"{adata.obs.outlier.value_counts()}")

    print("Outlier filtering using 3 x MAD for pct_counts_mt")
    adata.obs["mt_outlier"] = is_outlier(adata, "pct_counts_mt", 3) | (
            adata.obs["pct_counts_mt"] > 8
    )
    print(f"{adata.obs.mt_outlier.value_counts()}")

    print(f"Total number of cells: {adata.n_obs}")
    adata = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adata.n_obs}")

    ### run soupx
    print("## Running SoupX")
    soupx_path = os.path.join(qc_path, 'SoupX')
    if not os.path.exists(soupx_path):
        os.makedirs(soupx_path, exist_ok=True)

    adata = run_soupx(adata, soupx_path)

    # filter out genes that are not detected in at least 20 cells
    print("## Filter out genes that are not detected in at least 20 cells")
    print(f"Total number of genes: {adata.n_vars}")
    sc.pp.filter_genes(adata, min_cells=3)
    print(f"Number of genes after cell filter (> 20 cells): {adata.n_vars}")

    # filter cells with less than 100 genes expressed
    print("## Filter cells with less than 200 genes expressed")
    print(f"Total number of cells: {adata.n_obs}")
    sc.pp.filter_cells(adata, min_genes=200)
    print(f"Number of cells after filtering for min genes > 200: {adata.n_obs}")

    ## Doublet detection
    print("## Running Doublet detection")
    # run a doublet detection algorithm.Identifying doublets is crucial as they can lead to misclassifications or distortions
    # in downstream analysis steps.
    # TRY RUNNING scDblFinder R packjage in near future
    sc.pp.scrublet(adata, batch_key="sample")

    print("# Output QC graph after filtering for MT genes content < 5")
    output_qc(
        adata,
        os.path.join(qc_path, '02-AfterFiltering.pdf')
    )

    print("# Write resulting adata object to h5ad file")
    adata.write(os.path.join(outpath, 'adata_quality_control.h5ad'))
    return adata

def run_normalization(adata, outpath):
    print("## Running Normalization")
    ## Normalization
    # Count depth scaling normalizes the data to a “size factor” such as
    # 1 - the median count depth in the dataset,
    # 2 - ten thousand (CP10k)
    # 3 - or one million (CPM, counts per million).
    # The size factor for count depth scaling can be controlled via target_sum in pp.normalize_total.
    # We are applying median count depth normalization with log1p transformation (AKA log1PF).

    # shifted log method
    print("## Running shifted log method")
    scales_counts = sc.pp.normalize_total(adata, target_sum=None, inplace=False)
    # log1p transform
    adata.layers["log1p_norm"] = sc.pp.log1p(scales_counts["X"], copy=True)

    # Analytic Pearson residuals method
    print("## Analytic Pearson residuals method")
    analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(adata, inplace=False)
    adata.layers["analytic_pearson_residuals"] = csr_matrix(analytic_pearson["X"])

    print("# Write resulting adata object to h5ad file")
    adata.write(os.path.join(outpath, 'adata_normalisation.h5ad'))
    return adata

def run_feature_selection(adata, outpath):
    print("## Running Feature selection")
    fs_path = os.path.join(outpath, '02-FeatureSelection')
    if not os.path.exists(fs_path):
        os.makedirs(fs_path, exist_ok=True)

    #### seem to be not working well ####
    command = (
            "Rscript --vanilla "
            + os.path.join(os.path.dirname(__file__), '../R/run_scry.R ')
            + " -o " + fs_path
            + " -i " + os.path.join(outpath, 'adata_normalisation.h5ad')
    )
    print(command)
    return_code = subprocess.call(command, shell=True)

    if return_code == 0:
        print("Feature selection R script executed successfully.")
    else:
        print("Feature selection R script failed with return code", return_code)

    print("Loading Feature selection in adata")
    mask_df = pd.read_csv(os.path.join(fs_path, 'mask.tsv'), header=None, index_col=False, names=["highly_deviant"])
    mask = mask_df.iloc[:, 0]
    adata.var["highly_deviant"] = [m for m in mask]

    bd_df = pd.read_csv(os.path.join(fs_path, 'binomial_deviance.tsv'), header=None, index_col=False, names=["binomial_deviance"])
    binomial_deviance = bd_df.iloc[:, 0]
    adata.var["binomial_deviance"] = [b for b in binomial_deviance]

    # calculate dispersion field
    sc.pp.highly_variable_genes(adata, layer="log1p_norm", n_top_genes=2000)
    sc.experimental.pp.highly_variable_genes(
        adata, flavor="pearson_residuals", n_top_genes=2000
    )

    output_feature_selection_qc(adata, fs_path)

    print("# Write resulting adata object to h5ad file")
    adata.write(os.path.join(outpath, 'adata_feature_selection.h5ad'))

    return adata

def run_dim_reduction(adata, outpath):
    print("## Running Dimensionality Reduction")
    print(f"Number of features BEFORE dim reduction: {np.max(adata.X)}")

    dr_path = os.path.join(outpath, '03-DimReduction')
    if not os.path.exists(dr_path):
        os.makedirs(dr_path, exist_ok=True)

    adata.X = adata.layers["analytic_pearson_residuals"]
    # adata.var["highly_variable"] = adata.var["highly_deviant"]
    sc.pp.pca(adata, svd_solver="arpack", use_highly_variable=True)

    sc.tl.tsne(adata, use_rep="X_pca")
    sc.pp.neighbors(adata, n_pcs=30)
    sc.tl.umap(adata)

    output_dim_reduction_qc(adata, dr_path)

    print("# Write resulting adata object to h5ad file")
    adata.write(os.path.join(outpath, "adata_dimensionality_reduction.h5ad"))

    print(f"Number of features AFTER dim reduction: {np.max(adata.X)}")
    return adata

def output_dim_reduction_qc(adata, dr_path):
    print("Outputting Dim reduction plots")
    with PdfPages(os.path.join(dr_path, "DimReductionQC.pdf")) as pdf:

        sc.pl.pca(
            adata,
            color=["sample", "sample", "total_counts", "total_counts", "pct_counts_mt", "pct_counts_mt", "doublet_score", "doublet_score"],
            dimensions=[(0, 1), (2, 3), (0, 1), (2, 3), (0, 1), (2, 3), (0, 1), (2, 3)],
            ncols=4,
            size=3,
            show=False
        )
        pdf.savefig(bbox_inches="tight")

        sc.pl.tsne(
            adata,
            color=["sample", "total_counts", "pct_counts_mt", "doublet_score", "predicted_doublet"],
            ncols=5,
            show=False
        )
        pdf.savefig(bbox_inches="tight")

        sc.pl.umap(
            adata,
            color=["sample", "total_counts", "pct_counts_mt", "doublet_score", "predicted_doublet"],
            ncols=5,
            show=False
        )
        pdf.savefig(bbox_inches="tight")

        plt.close()

def run_clustering(adata, outpath, resolutions):
    print("## Running Clustering")
    cl_path = os.path.join(outpath, '04-Clustering')
    if not os.path.exists(cl_path):
        os.makedirs(cl_path, exist_ok=True)

    resKeys = ["sample"]
    for res in resolutions:
        print(f"Running Leiden clustering with resolution {res}")
        sc.tl.leiden(
            adata, key_added=f"leiden_{res}", resolution=res, flavor="igraph", n_iterations=2
        )
        resKeys.append(f"leiden_{res}")

    sc.pl.umap(
        adata,
        color=resKeys,
        ncols=5,
        show=False
    )
    plt.savefig(
        os.path.join(cl_path, 'Clustering.Allres.pdf'),
        bbox_inches="tight"
    )
    plt.close()

    print("# Write resulting adata object to h5ad file")
    adata.write(os.path.join(outpath, "adata_clustering.h5ad"))

    return adata

def run_celltypist_annotation(adata, ann_path, resolutions):
    print("## Running CellTypist Automated annotation")

    adata_celltypist = adata.copy()  # make a copy of our adata
    adata_celltypist.X = adata.layers["soupX_counts"]  # set adata.X to raw counts

    print("Normalise counts per cell to 10,000")
    sc.pp.normalize_per_cell(
        adata_celltypist, counts_per_cell_after=10 ** 4
    )  # normalize to 10,000 counts per cell
    sc.pp.log1p(adata_celltypist)  # log-transform
    # make .X dense instead of sparse, for compatibility with celltypist:
    adata_celltypist.X = adata_celltypist.X.toarray()

    print("Download celltypist models if not downloaded already: Immune_All_Low.pkl and Immune_All_High.pkl")
    models.download_models(
        force_update=True, model=["Immune_All_Low.pkl", "Immune_All_High.pkl"]
    )

    print("Load models")
    model_low = models.Model.load(model="Immune_All_Low.pkl")
    model_low.convert()
    model_low.write(os.path.join(ann_path, 'Immune_All_Low.mouse.pkl'))
    model_high = models.Model.load(model="Immune_All_High.pkl")
    model_high.convert()
    model_high.write(os.path.join(ann_path, 'Immune_All_High.mouse.pkl'))

    all_keys = []
    for res in resolutions:
        print(f"Running celltypist annotation with resolution {res}")
        print("Annotating using coarse model")
        cl_coarse = f"celltypist_cell_label_coarse_{res}"
        cs_coarse = f"celltypist_conf_score_coarse_{res}"
        predictions_high = celltypist.annotate(
            adata_celltypist, model=os.path.join(ann_path, 'Immune_All_High.mouse.pkl'), majority_voting=True, over_clustering=f"leiden_{res}"
        )
        predictions_high_adata = predictions_high.to_adata()

        adata.obs[cl_coarse] = predictions_high_adata.obs.loc[
            adata.obs.index, "majority_voting"
        ]
        adata.obs[cs_coarse] = predictions_high_adata.obs.loc[
            adata.obs.index, "conf_score"
        ]
        all_keys.append(cl_coarse)
        all_keys.append(cs_coarse)

        print("Annotating using fine model")
        cl_fine = f"celltypist_cell_label_fine_{res}"
        cs_fine = f"celltypist_conf_score_fine_{res}"
        predictions_low = celltypist.annotate(
            adata_celltypist, model=os.path.join(ann_path, 'Immune_All_Low.mouse.pkl'), majority_voting=True, over_clustering=f"leiden_{res}"
        )
        predictions_low_adata = predictions_low.to_adata()
        adata.obs[f"celltypist_cell_label_fine_{res}"] = predictions_low_adata.obs.loc[
            adata.obs.index, "majority_voting"
        ]
        adata.obs[f"celltypist_conf_score_fine_{res}"] = predictions_low_adata.obs.loc[
            adata.obs.index, "conf_score"
        ]
        all_keys.append(cl_fine)
        all_keys.append(cs_fine)


    output_annotation_qc(
        adata,
        os.path.join(ann_path, "CellTypist_Annotation.pdf"),
        all_keys
    )

    return adata


def output_annotation_qc(adata, o_path, all_keys):
    print("Outputting annotation plots")
    with PdfPages(o_path) as pdf:
        sc.pl.umap(
            adata,
            color=['sample'],
            frameon=False,
            sort_order=False,
            show=False
        )
        pdf.savefig(bbox_inches="tight")

        sc.pl.umap(
            adata,
            color=all_keys,
            frameon=False,
            sort_order=False,
            wspace=1,
            ncols=2,
            show=False
        )
        pdf.savefig(bbox_inches="tight")

        plt.close()

def run_clustifyr_annotation(adata, ann_path, resolutions):
    print("## Running Clustifyr annotation")

    print(f"Preparing clustifyr input files")
    adata_clustifyr = adata.copy()
    adata_clustifyr.X = adata.layers["soupX_counts"]
    data = adata_clustifyr.X.T
    sparse.save_npz(os.path.join(ann_path, 'data_matrix.npz'), data)
    pd.Series(adata_clustifyr.obs_names.to_series()).to_csv(
        os.path.join(ann_path, 'data_barcodes.tsv'),
        sep="\t", index=False, header=False
    )
    pd.Series(adata_clustifyr.var_names.to_series()).to_csv(
        os.path.join(ann_path, 'data_features.tsv'),
        sep="\t", index=False, header=False
    )

    highly_variant = adata_clustifyr.var["highly_variable"].loc[
        adata_clustifyr.var["highly_variable"] == True].index.tolist()
    with open(os.path.join(ann_path, 'highly_variant.tsv'), 'w', newline='\n') as f_output:
        tsv_output = csv.writer(f_output, delimiter='\n')
        tsv_output.writerow(highly_variant)

    all_keys = []
    for res in resolutions:
        print(f"Running clustifyr annotation with resolution {res}")
        clustifyr_groups = adata_clustifyr.obs[f"leiden_{res}"]
        clustifyr_groups.to_csv(os.path.join(ann_path, f'clusters_{res}.tsv'), sep="\t", index=True, header=False)

        print("Running Clustifyr in R")
        command = (
                "Rscript --vanilla "
                + os.path.join(os.path.dirname(__file__), '../R/run_clustifyr.R ')
                + " -o " + os.path.join(ann_path, f'clustifyr_leiden_{res}.tsv')
                + " -f " + os.path.join(ann_path, 'data_matrix.npz')
                + " -c " + os.path.join(ann_path, f'clusters_{res}.tsv')
                + " -g " + os.path.join(ann_path, 'data_features.tsv')
                + " -b " + os.path.join(ann_path, 'data_barcodes.tsv')
                + " -v " + os.path.join(ann_path, 'highly_variant.tsv')
                + " -x " + sys.executable
        )
        print(command)
        return_code = subprocess.call(command, shell=True)

        if return_code == 0:
            print("Clustifyr executed successfully.")
        else:
            print("Clustifyr failed with return code", return_code)

        print("Loading Clustifyr result matrix")
        clustifyr_matrix = pd.read_csv(os.path.join(ann_path, f'clustifyr_leiden_{res}.tsv'), sep="\t", header=None)

        celltype = clustifyr_matrix.iloc[:, 2]
        adata.obs[f"clustifyr_cell_label_{res}"] = [c for c in celltype]
        all_keys.append(f"clustifyr_cell_label_{res}")

        score = clustifyr_matrix.iloc[:, 3]
        adata.obs[f"clustifyr_conf_score_{res}"] = [c for c in score]
        all_keys.append(f"clustifyr_conf_score_{res}")

    output_annotation_qc(
        adata,
        os.path.join(ann_path, "Clustifyr_Annotation.pdf"),
        all_keys
    )

    return adata

# def run_TICAtlas_annotation(adata, ann_path, resolutions):
#     print("## Running TICAtlas annotation")
#     adata_to_map = adata.copy()
#     for layer in list(adata_to_map.layers.keys()):
#         if layer != "soupX_counts":
#             del adata_to_map.layers[layer]
#     adata_to_map.X = adata_to_map.layers["soupX_counts"]
#
#     # download from https://zenodo.org/records/5186413#.YRqbJC1h2v6
#     # ref: https://pmc.ncbi.nlm.nih.gov/articles/PMC8494216/
#     reference_model_features = pd.read_csv(
#         "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/TICAtlas_integrated_matrix.csv", index_col=0
#     )

def run_annotation(adata, outpath, resolutions, markers):
    ann_path = os.path.join(outpath, '05-Annotation')
    if not os.path.exists(ann_path):
        os.makedirs(ann_path, exist_ok=True)

    adata = run_celltypist_annotation(adata, ann_path, resolutions)
    adata = run_clustifyr_annotation(adata, ann_path, resolutions)
    # adata = run_TICAtlas_annotation(adata, ann_path, resolutions)

    for res in resolutions:
        print(f"Outputting markers annotation plots for {res}")

        with PdfPages(os.path.join(ann_path, f"Markers_Annotation_{res}.pdf"),) as pdf:

            adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X
            sc.pl.dotplot(
                adata,
                marker_genes,
                groupby=f"leiden_{res}",
                show=False
            )
            pdf.savefig(bbox_inches="tight")

            sc.pl.stacked_violin(
                adata,
                marker_genes,
                groupby=f"leiden_{res}",
                cmap="RdBu_r",
                show=False
            )
            pdf.savefig(bbox_inches="tight")

            # sc.pl.heatmap(
            #     adata,
            #     marker_genes,
            #     groupby=f"leiden_{res}",
            #     layer="scaled",
            #     dendrogram=True,
            #     show=False
            # )
            # pdf.savefig(bbox_inches="tight")

            sc.pl.matrixplot(
                adata,
                marker_genes,
                groupby=f"leiden_{res}",
                colorbar_title="mean z-score",
                layer="scaled",
                vmin=-2,
                vmax=2,
                cmap="RdBu_r",
                show=False,
            )
            pdf.savefig(bbox_inches="tight")

            plt.close()

    print("# Write resulting adata object to h5ad file")
    adata.write(os.path.join(outpath, "adata_annotation.h5ad"))

    return adata

def aggregate_and_filter(
    adata,
    cell_identity,
    donor_key="sample",
    cell_identity_key="cell_type",
    obs_to_keep=[],  # which additional metadata to keep, e.g. gender, age, etc.
    replicates_per_patient=3,
):
    NUM_OF_CELL_PER_DONOR = 30
    # subset adata to the given cell identity
    adata_cell_pop = adata[adata.obs[cell_identity_key] == cell_identity].copy()
    # check which donors to keep according to the number of cells specified with NUM_OF_CELL_PER_DONOR
    size_by_donor = adata_cell_pop.obs.groupby([donor_key]).size()
    donors_to_drop = [
        donor
        for donor in size_by_donor.index
        if size_by_donor[donor] <= NUM_OF_CELL_PER_DONOR
    ]
    if len(donors_to_drop) > 0:
        print("Dropping the following samples:")
        print(donors_to_drop)
    df = pd.DataFrame(columns=[*adata_cell_pop.var_names, *obs_to_keep])

    adata_cell_pop.obs[donor_key] = adata_cell_pop.obs[donor_key].astype("category")
    for i, donor in enumerate(donors := adata_cell_pop.obs[donor_key].cat.categories):
        print(f"\tProcessing donor {i+1} out of {len(donors)}...", end="\r")
        if donor not in donors_to_drop:
            adata_donor = adata_cell_pop[adata_cell_pop.obs[donor_key] == donor]
            # create replicates for each donor
            indices = list(adata_donor.obs_names)
            random.shuffle(indices)
            indices = np.array_split(np.array(indices), replicates_per_patient)
            for i, rep_idx in enumerate(indices):
                adata_replicate = adata_donor[rep_idx]
                # specify how to aggregate: sum gene expression for each gene for each donor and also keep the condition information
                agg_dict = {gene: "sum" for gene in adata_replicate.var_names}
                for obs in obs_to_keep:
                    agg_dict[obs] = "first"
                # create a df with all genes, donor and condition info
                df_donor = pd.DataFrame(adata_replicate.X.A)
                df_donor.index = adata_replicate.obs_names
                df_donor.columns = adata_replicate.var_names
                df_donor = df_donor.join(adata_replicate.obs[obs_to_keep])
                # aggregate
                df_donor = df_donor.groupby(donor_key).agg(agg_dict)
                df_donor[donor_key] = donor
                df.loc[f"donor_{donor}_{i}"] = df_donor.loc[donor]
    print("\n")
    # create AnnData object from the df
    adata_cell_pop = sc.AnnData(
        df[adata_cell_pop.var_names], obs=df.drop(columns=adata_cell_pop.var_names)
    )
    return adata_cell_pop

def run_dge(adata, outpath, resolutions):
    print("## Running differential gene expression analysis")
    dge_path = os.path.join(outpath, '06-DifferentialGeneExpression')
    if not os.path.exists(dge_path):
        os.makedirs(dge_path, exist_ok=True)

    # x = np.max(adata.X)
    print(f"Preparing dge input files")
    adata_dge = adata.copy()
    adata_dge.layers["counts"] = adata_dge.layers["soupX_counts"]

    adata_dge.obs["clustifyr_cell_label_0.1"] = [ct.replace(" ", "_") for ct in adata_dge.obs["clustifyr_cell_label_0.1"]]
    adata_dge.obs["clustifyr_cell_label_0.1"] = [ct.replace("+", "") for ct in adata_dge.obs["clustifyr_cell_label_0.1"]]

    adata_dge.obs["sample"] = adata_dge.obs["sample"].astype("category")
    adata_dge.obs["clustifyr_cell_label_0.1"] = adata_dge.obs["clustifyr_cell_label_0.1"].astype("category")

    adata_dge.X = adata_dge.layers["counts"].copy()
    obs_to_keep = ["cell_type", "sample"]

    cell_type = adata_dge.obs["clustifyr_cell_label_0.1"].cat.categories[0]
    print(
        f'Processing {cell_type} (1 out of {len(adata.obs["clustifyr_cell_label_0.1"].cat.categories)})...'
    )
    adata_pb = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)
    for i, cell_type in enumerate(adata.obs["cell_type"].cat.categories[1:]):
        print(
            f'Processing {cell_type} ({i + 2} out of {len(adata.obs["cell_type"].cat.categories)})...'
        )
        adata_cell_type = aggregate_and_filter(adata, cell_type, obs_to_keep=obs_to_keep)
        adata_pb = adata_pb.concatenate(adata_cell_type)

    print("done dge")


def main(samples, outpath, resolutions, markers):
    print("Starting")
    # adata = read_cellranger(samples)
    #
    # # create out dir if not exists
    # if not os.path.exists(outpath):
    #     os.makedirs(outpath, exist_ok=True)
    #
    # adata = run_qc(adata, outpath)
    #
    # # adata = sc.read_h5ad(os.path.join(outpath, "adata_quality_control.h5ad"))
    # adata = run_normalization(adata, outpath)
    #
    # # adata = sc.read_h5ad(os.path.join(outpath, "adata_normalisation.h5ad"))
    # adata = run_feature_selection(adata, outpath)
    #
    # # adata = sc.read_h5ad(os.path.join(outpath, "adata_feature_selection.h5ad"))
    # adata = run_dim_reduction(adata, outpath)
    #
    # # adata = sc.read_h5ad(os.path.join(outpath, "adata_dimensionality_reduction.h5ad"))
    # adata = run_clustering(adata, outpath, resolutions)

    # # adata = sc.read_h5ad(os.path.join(outpath, "adata_clustering.h5ad"))
    # adata = run_annotation(adata, outpath, resolutions, markers)

    adata = sc.read_h5ad(os.path.join(outpath, "adata_annotation.h5ad"))
    adata = run_dge(adata, outpath, resolutions)

    print("Done")

samples = {
    "FL_T": "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/cellranger/FL_T/filtered_feature_bc_matrix.h5",
    "V_T": "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/cellranger/V_T/filtered_feature_bc_matrix.h5",
}
samples_raw = {
    "FL_T": "/home/jflucier/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/cellranger/FL_T/raw_feature_bc_matrix.h5",
    "V_T": "/home/jflucier/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/cellranger/V_T/raw_feature_bc_matrix.h5",
}

# leiden resolutions to process
resolutions = [0.10, 0.25, 0.50, 1.00]

## Marker gene set
marker_genes = {
    "CD45": ["Ptprc"],
    "CD8+ T": ["Cd8a", "Cd8b1", "Gzmk", "Nkg7"],
    "CD4+ T": ["Cd4", "Ctla4", "Icos"],
    "CD3+ T": ["Cd3d", "Cd3e", "Cd3g", "Lck"],
    "B cells": ["Cd79a", "Cd79b", "Ighm"], #"Cd19",
    "NK cells": ["Gzma", "Klra4", "Nkg7", "Prf1"],
    "Macrophages": ["Adgre1", "Cd14", "Fcgr3", "Mafb"],
    "Monocytes": ["Ly6c2", "Spn", "Ccl2", "Ccl7"],
    "cDCs": ["Xcr1", "Flt3", "Ccr7", "Cd86"],
    "pDCs": ["Siglech", "Clec10a", "Clec12a", "Lair1"],
    "Neutrophils": ["Csf3r", "Cxcl3"] #, "S100a8", "S100a9"
}

main(samples, outpath, resolutions, marker_genes)




# # Let us compute the neighborhood graph of cells using the PCA representation of the data matrix.
# # sc.pp.neighbors(adata)
# # or used bknn
# print("# BKNN")
# sc.external.pp.bbknn(adata, batch_key="sample")



# ## Manual cell-type annotation
# for res in resKeys:
#
#     if res == "sample":
#         continue
#
#     # output count column by res
#     counts_t = getattr(adata.obs, res).value_counts()
#     counts_t.rename(f"total_count", inplace=True)
#     list_of_series = [counts_t]
#     cols=['total_counts']
#     for sample_id, filename in samples.items():
#         adata_s = adata[adata.obs[f"sample"] == sample_id]
#         counts = getattr(adata_s.obs, res).value_counts()
#         counts.rename(f"{sample_id}_count", inplace=True)
#         list_of_series.append(counts)
#         # cols.append(f"{sample_id}_count")
#
#     # df = pd.DataFrame(list_of_series, columns=cols, index=counts_t.index)
#     df = pd.concat(list_of_series, axis=1)
#     df.to_csv(
#         f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.counts.tsv",
#         header=True,
#         index=True,
#         sep="\t"
#     )
#
#
#     print(f"## Output UMAP for for {res}")
#     adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X
#     sc.pl.umap(
#         adata,
#         color=["sample", f"{res}"],
#         show=False
#     )
#     plt.savefig(
#         f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.umap.pdf",
#         bbox_inches="tight"
#     )
#     plt.close()
#
#     print(f"## Running Differentially-expressed Genes as Markers for {res}")
#


#     # Obtain cluster-specific differentially expressed genes
#     sc.tl.rank_genes_groups(adata, groupby=f"{res}", method="wilcoxon")
#     df = sc.get.rank_genes_groups_df(adata, group=None,  pval_cutoff=0.01)
#     df.to_csv(
#         f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/DGE/FL_T_vs_V_T.{res}.dge.tsv",
#         header=True,
#         index=False,
#         sep="\t"
#     )
#
#     # We can then visualize the top 5 differentially-expressed genes on a dotplot.
#     for n in [5, 10, 25]:
#         with PdfPages(f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.intercluster.top{n}.pdf") as pdf:
#             sc.pl.rank_genes_groups_dotplot(
#                 adata,
#                 groupby=f"{res}",
#                 # standard_scale="var",
#                 n_genes=n,
#                 values_to_plot="logfoldchanges",
#                 vmax=7,
#                 vmin=-7,
#                 cmap="bwr",
#                 return_fig=True, show=False
#             )
#
#             sc.pl.rank_genes_groups_matrixplot(
#                 adata,
#                 groupby=f"{res}",
#                 n_genes=n,
#                 use_raw=False,
#                 vmin=-3, vmax=3,
#                 cmap="bwr",
#                 layer="scaled",
#                 colorbar_title="mean z-score",
#                 return_fig=True, show=False
#             )
#
#             sc.pl.rank_genes_groups(adata, n_genes=n, sharey=False, show=False,  return_fig=True)
#
#             dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group=None, key="wilcoxon").head(n)["names"]
#             sc.pl.umap(
#                 adata,
#                 color=[f"{res}", "sample", *dc_cluster_genes],
#                 frameon=False,
#                 ncols=3,
#                 show=False,
#                 return_fig=True
#             )
#
#         sc.pl.rank_genes_groups_tracksplot(
#             adata,
#             groupby=f"{res}",
#             n_genes=n,
#             show=False, return_fig=True
#         )
#
#         sc.pl.rank_genes_groups_violin(
#             adata,
#             n_genes=n,
#             jitter=False,
#             show=False
#         )



