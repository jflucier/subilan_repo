# Core scverse libraries
import scanpy as sc
import anndata as ad
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd

sc.settings.set_figure_params(dpi=50, facecolor="white")

samples = {
    "FL_T": "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/cellranger/FL_T/filtered_feature_bc_matrix.h5",
    "V_T": "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/cellranger/V_T/filtered_feature_bc_matrix.h5",
}
adatas = {}

print("Reading cellranger count matrix")
for sample_id, filename in samples.items():
    # path = EXAMPLE_DATA.fetch(filename)
    sample_adata = sc.read_10x_h5(filename)
    sample_adata.var_names_make_unique()
    adatas[sample_id] = sample_adata

adata = ad.concat(adatas, label="sample")
adata.obs_names_make_unique()
# print(adata.obs["sample"].value_counts())

## Quality Control
print("## Running QC")
# mitochondrial genes, "MT-" for human, "Mt-" for mouse
adata.var["mt"] = adata.var_names.str.startswith("mt-")
# ribosomal genes
adata.var["ribo"] = adata.var_names.str.startswith(("Rps", "Rpl"))

sc.pp.calculate_qc_metrics(
    adata, qc_vars=["mt", "ribo"], inplace=True, log1p=True
)

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
plt.savefig("/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/QC/FL_T_vs_V_T.beforeqc.pdf", bbox_inches="tight")

# scatter plot colored by pct_counts_mt
sc.pl.scatter(
    adata,
    "total_counts",
    "n_genes_by_counts",
    color="pct_counts_mt",
    show=False
)
plt.savefig("/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/QC/FL_T_vs_V_T.beforeqc.MT.pdf", bbox_inches="tight")
plt.close()
sc.pl.scatter(
    adata,
    "total_counts",
    "n_genes_by_counts",
    color="pct_counts_ribo",
    show=False
)
plt.savefig("/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/QC/FL_T_vs_V_T.beforeqc.RIBO.pdf", bbox_inches="tight")
plt.close()

# filter cells with less than 100 genes expressed
sc.pp.filter_cells(adata, min_genes=100)
# genes that are detected in less than 3 cells
sc.pp.filter_genes(adata, min_cells=3)

# filter for percent mito
# High proportions are indicative of poor-quality cells (Islam et al. 2014; Ilicic et al. 2016),
# possibly because of loss of cytoplasmic RNA from perforated cells. The reasoning is that mitochondria are larger
# than individual transcript molecules and less likely to escape through tears in the cell membrane
print("Number cells before filtering for MT genes content < 5: %d"%adata.n_obs)
adata = adata[adata.obs['pct_counts_mt'] < 5, :]
print("Remaining cells: %d"%adata.n_obs)
# filter for percent ribo > 0.05
# print("Number cells before filtering for ribosomal genes content > 5: %d"%adata.n_obs)
# adata = adata[adata.obs['pct_counts_ribo'] > 5, :]
# print("Remaining cells: %d"%adata.n_obs)

sc.pl.violin(
    adata,
    ["n_genes_by_counts", "total_counts", "pct_counts_mt", "pct_counts_ribo"],
    jitter=0.4,
    multi_panel=True,
    show=False
)
plt.savefig("/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/QC/FL_T_vs_V_T.afterqc.pdf", bbox_inches="tight")
plt.close()

# scatter plot colored by pct_counts_mt
sc.pl.scatter(
    adata,
    "total_counts",
    "n_genes_by_counts",
    color="pct_counts_mt",
    show=False
)
plt.savefig("/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/QC/FL_T_vs_V_T.afterqc.MT.pdf", bbox_inches="tight")
plt.close()
sc.pl.scatter(
    adata,
    "total_counts",
    "n_genes_by_counts",
    color="pct_counts_ribo",
    show=False
)
plt.savefig("/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/QC/FL_T_vs_V_T.afterqc.RIBO.pdf", bbox_inches="tight")
plt.close()

## Doublet detection
print("## Running Doublet detection")
# run a doublet detection algorithm.Identifying doublets is crucial as they can lead to misclassifications or distortions
# in downstream analysis steps.
# adds doublet_score and predicted_doublet to .obs
# can now either filter directly on predicted_doublet or use the doublet_score later during clustering
# to filter clusters with high doublet scores.
sc.pp.scrublet(adata, batch_key="sample")

## Normalization
# Count depth scaling normalizes the data to a “size factor” such as
# 1 - the median count depth in the dataset,
# 2 - ten thousand (CP10k)
# 3 - or one million (CPM, counts per million).
# The size factor for count depth scaling can be controlled via target_sum in pp.normalize_total.
# We are applying median count depth normalization with log1p transformation (AKA log1PF).

# Saving count data
adata.layers["counts"] = adata.X.copy()
# Normalizing to median total counts
# If adding parameter `target_sum=1e6`, this is CPM normalization.
# If adding parameter `target_sum=1e4`, this is CP10k.
sc.pp.normalize_total(adata, target_sum=1e4)
# Logarithmize the data
sc.pp.log1p(adata)

## Feature selection
print("## Running Feature selection")
# The scanpy function pp.highly_variable_genes annotates highly variable genes by reproducing the implementations
# of Seurat [Satija et al., 2015],
# Cell Ranger [Zheng et al., 2017],
# and Seurat v3 [Stuart et al., 2019] depending on the chosen flavor.

# param flavor: Literal["seurat", "cell_ranger", "seurat_v3", "seurat_v3_paper"]. Default = "seurat"
# `'seurat_v3'`/`'seurat_v3_paper'` requires `scikit-misc` package.
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="sample")
sc.pl.highly_variable_genes(adata, show=False)
plt.savefig(
    "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/QC/FL_T_vs_V_T.highly_variable_genes.pdf",
    bbox_inches="tight"
)
plt.close()

## Dimensionality Reduction
print("## Running Dimensionality Reduction")
# Reduce the dimensionality of the data by running principal component analysis (PCA)
sc.tl.pca(adata)
# Let us inspect the contribution of single PCs to the total variance in the data.
# This gives us information about how many PCs we should consider in order to compute the
# neighborhood relations of cells, e.g. used in the clustering function leiden() or tsne().
# In our experience, there does not seem to be signifigant downside to overestimating the numer of principal components.
sc.pl.pca_variance_ratio(adata, n_pcs=50, log=True, show=False)
plt.savefig(
    "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/QC/FL_T_vs_V_T.dim_reduction.pdf",
    bbox_inches="tight"
)
plt.close()

# You can also plot the principal components to see if there are any potentially undesired features (e.g. batch, QC metrics)
# driving signifigant variation in this dataset.
# In this case, there isn’t anything too alarming, but it’s a good idea to explore this.
sc.pl.pca(
    adata,
    color=["sample", "sample", "pct_counts_mt", "pct_counts_mt"],
    dimensions=[(0, 1), (2, 3), (0, 1), (2, 3)],
    ncols=2,
    size=2,
    show=False
)
plt.savefig(
    "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/QC/FL_T_vs_V_T.pca.pdf",
    bbox_inches="tight"
)
plt.close()

## Nearest neighbor graph constuction and visualization
print("## Running Nearest neighbor graph constuction and visualization")
# Let us compute the neighborhood graph of cells using the PCA representation of the data matrix.
sc.pp.neighbors(adata)
# This graph can then be embedded in two dimensions for visualiztion with UMAP
sc.tl.umap(adata)
# We can now visualize the UMAP according to the sample.
sc.pl.umap(
    adata,
    color="sample",
    # Setting a smaller point size to get prevent overlap
    size=2,
    show=False
)
plt.savefig(
    "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.samples.pdf",
    bbox_inches="tight"
)
plt.close()

### If you inspect batch effects in your UMAP it can be beneficial to integrate across samples
### and perform batch correction/integration. We recommend checking out scanorama and scvi-tools for batch integration.

## Clustering
print("## Running Clustering")
# we recommend the Leiden graph-clustering method (community detection based on optimizing modularity)
# Using the igraph implementation and a fixed number of iterations can be significantly faster, especially for larger datasets
# flavor: Literal["leidenalg", "igraph"] = "leidenalg"
# n_iterations: -1 has the algorithm run until it reaches its optimal clustering. 2 is faster and the default for underlying packages.
sc.tl.leiden(adata, flavor="igraph", n_iterations=2)
sc.pl.umap(adata, color=["leiden"], show=False)
plt.savefig(
    "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.clusters.pdf",
    bbox_inches="tight"
)
plt.close()

## Re-assess quality control and cell filtering
# visualizing different QC metrics using UMAP.
sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
    show=False
)
plt.savefig(
    "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.doublets.before.pdf",
    bbox_inches="tight"
)
plt.close()

print("Number cells before filtering for doublets: %d"%adata.n_obs)
adata = adata[adata.obs['predicted_doublet'] == False, :]
print("Remaining cells: %d"%adata.n_obs)

sc.pl.umap(
    adata,
    color=["leiden", "predicted_doublet", "doublet_score"],
    # increase horizontal space between panels
    wspace=0.5,
    size=3,
    show=False
)
plt.savefig(
    "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.doublets.after.pdf",
    bbox_inches="tight"
)
plt.close()

sc.pl.umap(
    adata,
    color=["leiden", "log1p_total_counts", "pct_counts_mt", "pct_counts_ribo", "log1p_n_genes_by_counts"],
    wspace=0.5,
    ncols=2,
    show=False
)
plt.savefig(
    "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.qc_metrics.pdf",
    bbox_inches="tight"
)
plt.close()

## Manual cell-type annotation
print("## Running Manual cell-type annotation")
# let’s generate a set of clustering solutions which we can then use to annotate our cell types.
# resolution: A parameter value controlling the coarseness of the clustering. Higher values lead to more clusters.
# resArr = [0.1, 0.2, 0.3, 0.5, 0.75, 1.0]
resArr = [0.2]
resKeys = ["sample"]
for res in resArr:
    sc.tl.leiden(
        adata, key_added=f"leiden_{res}", resolution=res, flavor="igraph"
    )
    resKeys.append(f"leiden_{res}")

sc.pl.umap(
    adata,
    color=resKeys,
    show=False
)
plt.savefig(
    "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.allresolutions.umap.pdf",
    bbox_inches="tight"
)
plt.close()

## Marker gene set
marker_genes = {
    "CD8+ T": ["Cd8a", "Cd8b1"],
    "CD4+ T": ["Cd4", "Ctla4", "Icos"],
    "CD3+ ab- T": ["Cd3d", "Cd3e", "Cd3g", "Lck"],
    "B cells": ["Cd79a", "Cd79b"],
    "NK cells": ["Gzma", "Klra4", "Nkg7", "Prf1"],
    "Macrophages": ["Adgre1", "Cd14", "Fcgr3", "Mafb"],
    "Monocytes": ["Ly6c2", "Spn", "Ccl2", "Ccl7"],
    "cDCs": ["Xcr1", "Flt3", "Ccr7", "Cd86"],
    "pDCs": ["Siglech", "Clec10a", "Clec12a"],
    "Neutrophils": ["Csf3r", "Cxcl3"],
    "CD45": ["Ptprc"]
}

for res in resKeys:

    if res == "sample":
        continue

    # output count column by res
    counts_t = getattr(adata.obs, res).value_counts()
    counts_t.rename(f"total_count", inplace=True)
    list_of_series = [counts_t]
    cols=['total_counts']
    for sample_id, filename in samples.items():
        adata_s = adata[adata.obs[f"sample"] == sample_id]
        counts = getattr(adata_s.obs, res).value_counts()
        counts.rename(f"{sample_id}_count", inplace=True)
        list_of_series.append(counts)
        # cols.append(f"{sample_id}_count")

    # df = pd.DataFrame(list_of_series, columns=cols, index=counts_t.index)
    df = pd.concat(list_of_series, axis=1)
    df.to_csv(
        f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.counts.tsv",
        header=True,
        index=True,
        sep="\t"
    )


    print(f"## Output UMAP for for {res}")
    adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X
    sc.pl.umap(
        adata,
        color=["sample", f"{res}"],
        show=False
    )
    plt.savefig(
        f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.umap.pdf",
        bbox_inches="tight"
    )
    plt.close()

    print(f"## Running Differentially-expressed Genes as Markers for {res}")

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(15, 15))

    sc.pl.dotplot(
        adata,
        marker_genes,
        groupby=f"{res}",
        # colorbar_title="mean z-score",
        # layer="scaled",
        # vmin=-2,
        # vmax=2,
        cmap="RdBu_r",
        ax = ax1,
        show=False
    )
    # plt.savefig(
    #     f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.markers.pdf",
    #     bbox_inches="tight"
    # )
    # plt.close()

    sc.pl.stacked_violin(
        adata,
        marker_genes,
        groupby=f"{res}",
        cmap="RdBu_r",
        ax=ax2,
        show=False
    )
    # plt.savefig(
    #     f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.markers.violin.pdf",
    #     bbox_inches="tight"
    # )
    # plt.close()

    # sc.pl.heatmap(
    #     adata,
    #     marker_genes,
    #     groupby=f"{res}",
    #     cmap="RdBu_r",
    #     dendrogram=True,
    #     ax=ax3,
    #     show=False
    # )

    # sc.pl.matrixplot(
    #     adata,
    #     marker_genes,
    #     groupby=f"{res}",
    #     colorbar_title="mean z-score",
    #     layer="scaled",
    #     vmin=-2,
    #     vmax=2,
    #     cmap="RdBu_r",
    #     ax=ax3,
    #     show=False,
    # )

    plt.savefig(
        f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.markers.pdf",
        bbox_inches="tight"
    )
    plt.close()

    # adata.layers["scaled"] = sc.pp.scale(adata, copy=True).X
    # sc.pl.matrixplot(
    #     adata,
    #     marker_genes,
    #     groupby=f"{res}",
    #     dendrogram=True,
    #     colorbar_title="mean z-score",
    #     layer="scaled",
    #     vmin=-2,
    #     vmax=2,
    #     cmap="RdBu_r",
    #     show=False,
    # )
    # plt.savefig(
    #     f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.markers.matrix.pdf",
    #     bbox_inches="tight"
    # )
    # plt.close()

    # Obtain cluster-specific differentially expressed genes
    sc.tl.rank_genes_groups(adata, groupby=f"{res}", method="wilcoxon", key_added = "wilcoxon")
    df = sc.get.rank_genes_groups_df(adata, group=None, key=f"wilcoxon", pval_cutoff=0.01)
    df.to_csv(
        f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/DGE/FL_T_vs_V_T.{res}.dge.tsv",
        header=True,
        index=False,
        sep="\t"
    )

    sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False, show=False, key="wilcoxon")
    plt.savefig(
        f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/DGE/FL_T_vs_V_T.{res}.dge.top25.pdf",
        bbox_inches="tight"
    )
    plt.close()

    # We can then visualize the top 5 differentially-expressed genes on a dotplot.
    for n in [5, 10, 25]:
        with PdfPages(f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.intercluster.top{n}.pdf") as pdf:
            fig = sc.pl.rank_genes_groups_dotplot(
                adata, groupby=f"{res}", standard_scale="var", n_genes=n, key="wilcoxon", return_fig=True
            )
            pdf.savefig(fig)
            # plt.savefig(
            #     f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.intercluster.top{n}.genes.pdf",
            #     bbox_inches="tight"
            # )
            # plt.close()

            dc_cluster_genes = sc.get.rank_genes_groups_df(adata, group=None, key="wilcoxon").head(n)["names"]
            sc.pl.umap(
                adata,
                color=[*dc_cluster_genes, f"{res}", "sample"],
                frameon=False,
                ncols=3,
                show=False
            )
            pdf.savefig(fig)
            # plt.savefig(
            #     f"/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/UMAP/FL_T_vs_V_T.{res}.intercluster.top{n}.pdf",
            #     bbox_inches="tight"
            # )
            # plt.close()



# adata.obs["cell_type_lvl1"] = adata.obs["leiden_res_0.20"].map(
#     {
#         "0": "CD8+ T",
#         "1": "CD8++ T",
#         "2": "CD8++ T",
#         "3": "CD4+ T",
#         "4": "CD4- CD8- T",
#     }
# )
# sc.pl.dotplot(adata, marker_genes, groupby="leiden_res_1.00", standard_scale="var", show=False)
# plt.savefig(
#     "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/FL_T_vs_V_T.qc.umapclusterbymarkers.celltype.pdf",
#     bbox_inches="tight"
# )
# plt.close()

## Differentially-expressed Genes as Markers
# print("## Running Differentially-expressed Genes as Markers")
# # Obtain cluster-specific differentially expressed genes
# sc.tl.rank_genes_groups(adata, groupby="leiden_res_1.00", method="wilcoxon")
#
# # We can then visualize the top 5 differentially-expressed genes on a dotplot.
# sc.pl.rank_genes_groups_dotplot(
#     adata, groupby="leiden_res_1.00", standard_scale="var", n_genes=10, show=False
# )
# plt.savefig(
#     "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/scanpy/FL_T_vs_V_T.qc.dgebygroup.top10.pdf",
#     bbox_inches="tight"
# )
# plt.close()
