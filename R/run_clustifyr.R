
# clustifyr
# https://github.com/rnabioco/clustifyr

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("clustifyr")
# BiocManager::install("clustifyrdatahub")
# BiocManager::install("SummarizedExperiment")

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("clustifyr"))
suppressPackageStartupMessages(library("clustifyrdatahub"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("reticulate"))
scipy_sparse <- import("scipy.sparse")
sc <- import("scanpy")
np <- import("numpy")

option_list <- list(
  make_option(c("-o", "--out"), type="character", default=NULL, help="output file path", metavar="character"),
  make_option(c("-f", "--filtered_npz"), type="character", default=NULL, help="Input files filtered data matrix npz file", metavar="character"),
  make_option(c("-g", "--features"), type="character", default=NULL, help="File data_features.tsv, list of genes", metavar="character"),
  make_option(c("-b", "--barcodes"), type="character", default=NULL, help="File data_barcodes.tsv, list of cellss", metavar="character"),
  make_option(c("-c", "--clusters"), type="character", default=NULL, help="Input clusters.tsv, list of cells<tab>cluster_id", metavar="character"),
  make_option(c("-v", "--highly_variant"), type="character", default=NULL, help="Input highly variable gene list", metavar="character"),
  make_option(c("-x", "--python_path"), type="character", default="/storage/Documents/service/biologie/venv/bin/python", help="Path to python exec", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

out_f <- opt$out
filtered_npz <- opt$filtered_npz
raw_npz <- opt$raw_npz
features <- opt$features
barcodes <- opt$barcodes
clusters <- opt$clusters
highly_variant <- opt$highly_variant
python <- opt$python_path

# out_f <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/05-Annotation/clustifyr_lenden_0.25.tsv"
# filtered_npz <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/05-Annotation/data_matrix.npz"
# features <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/05-Annotation/data_features.tsv"
# barcodes <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/05-Annotation/data_barcodes.tsv"
# clusters <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/05-Annotation/clusters.tsv"
# highly_variant <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/05-Annotation/highly_variant.tsv"
# python <- "/home/jflucier/Documents/service/biologie/venv/bin/python"

print(paste0("out_f=",out_f))
print(paste0("filtered_npz=",filtered_npz))
print(paste0("features=",features))
print(paste0("barcodes=",barcodes))
print(paste0("clusters=",clusters))
print(paste0("highly_variant=",highly_variant))
print(paste0("python=",python))

use_python(python)

out_dir <- dirname(out_f)
setwd(out_dir)

print("loading filtered data matrix")
data <- scipy_sparse$load_npz(filtered_npz)

print("loading gene list")
g_data <- read.csv(
  features,
  header = FALSE,
  sep="\t",
  na.strings = "",
  stringsAsFactors=FALSE
)

print("loading cells list")
b_data <- read.csv(
  barcodes,
  header = FALSE,
  sep="\t",
  na.strings = "",
  stringsAsFactors=FALSE
)

print("loading clusters")
cl <- read.csv(
  clusters,
  header = FALSE,
  col.names = c("cellid", "cluster"),
  row.names = 1,
  sep="\t",
  na.strings = "",
  stringsAsFactors=FALSE
)

print("loading highly variant genes")
hv <- read.csv(
  highly_variant,
  header = FALSE,
  sep="\t",
  na.strings = "",
  stringsAsFactors=FALSE
)


nv_cl <- as.data.frame(as.character(cl[["cluster"]]))
colnames(nv_cl) <- c("cluster")
rownames(data) <- g_data[["V1"]]
colnames(data) <- b_data[["V1"]]
hvg <- hv[["V1"]]
data <- as(data, "sparseMatrix")

print("retreving reference immgen")
ref <- ref_immgen()

print("Running clustifyr")
print("assign celltype using immgen reference")
res <- clustify(
  input = data,
  metadata = nv_cl,
  cluster_col = "cluster",
  ref_mat = ref,
  query_genes = hvg,
  per_cell = FALSE
)

x <- cor_to_call(res)
cl$type <- NA
cl$r <- NA
cl$type <- x$type[match(cl$cluster, x$cluster)]
cl$r <- x$r[match(cl$cluster, x$cluster)]


print(paste0("Outputting results in ",out_dir))
write.table(
  cl,
  file = out_f,
  sep = "\t",
  row.names = TRUE,
  col.names = FALSE
)

print("Done!")