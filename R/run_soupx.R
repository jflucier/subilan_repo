
# SoupX
# https://rawcdn.githack.com/constantAmateur/SoupX/204b602418df12e9fdb4b68775a8b486c6504fe4/inst/doc/pbmcTutorial.html

# # install.packages("SoupX")
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# install.packages("reticulate")

suppressPackageStartupMessages(library("Matrix"))
suppressPackageStartupMessages(library("SoupX"))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("reticulate"))
scipy_sparse <- import("scipy.sparse")

option_list <- list(
  make_option(c("-o", "--out_npz"), type="character", default=NULL, help="output file path", metavar="character"),
  make_option(c("-f", "--filtered_npz"), type="character", default=NULL, help="Input files filtered data matrix npz file", metavar="character"),
  make_option(c("-r", "--raw_npz"), type="character", default=NULL, help="Input files raw data matrix npz file", metavar="character"),
  make_option(c("-g", "--features"), type="character", default=NULL, help="File data_features.tsv, list of genes", metavar="character"),
  make_option(c("-b", "--barcodes"), type="character", default=NULL, help="File data_barcodes.tsv, list of cellss", metavar="character"),
  make_option(c("-c", "--clusters"), type="character", default=NULL, help="Input clusters.tsv, list of cells<tab>cluster_id", metavar="character"),
  make_option(c("-x", "--python_path"), type="character", default="/storage/Documents/service/biologie/venv/bin/python", help="Path to python exec", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

out_npz <- opt$out_npz
filtered_npz <- opt$filtered_npz
raw_npz <- opt$raw_npz
features <- opt$features
barcodes <- opt$barcodes
clusters <- opt$clusters
python <- opt$python_path

# out_npz <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/SoupX/soupx_matrix.npz"
# filtered_npz <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/SoupX/data_matrix.npz"
# raw_npz <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/SoupX/data_tod_matrix.npz"
# features <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/SoupX/data_features.tsv"
# barcodes <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/SoupX/data_barcodes.tsv"
# clusters <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/SoupX/clusters.tsv"
# python <- "/home/jflucier/Documents/service/biologie/venv/bin/python"

print(paste0("out_npz=",out_npz))
print(paste0("filtered_npz=",filtered_npz))
print(paste0("raw_npz=",raw_npz))
print(paste0("features=",features))
print(paste0("barcodes=",barcodes))
print(paste0("clusters=",clusters))
print(paste0("python=",python))


use_python(python)

out_dir <- dirname(out_npz)
setwd(out_dir)

print("loading filtered data matrix")
data <- scipy_sparse$load_npz(filtered_npz)
print("loading raw data matrix")
data_tod <- scipy_sparse$load_npz(raw_npz)

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
  row.names = 1,
  sep="\t",
  na.strings = "",
  stringsAsFactors=FALSE
)

nv_cl <- as.character(cl[["V2"]])
names(nv_cl) <- rownames(cl)
rownames(data) <- g_data[["V1"]]
colnames(data) <- b_data[["V1"]]
data <- as(data, "sparseMatrix")
data_tod <- as(data_tod, "sparseMatrix")

print("Generate SoupChannel Object for SoupX")
sc <- SoupChannel(data_tod, data, calcSoupProfile = FALSE)

print("Add extra meta data to the SoupChannel object")
soupProf <- data.frame(row.names = rownames(data), est = rowSums(data)/sum(data), counts = rowSums(data))
sc <- setSoupProfile(sc, soupProf)
print("Set cluster information in SoupChannel")
sc <- setClusters(sc, nv_cl)

print("Estimate contamination fraction")
png(paste0(out_dir,"/autoEstCont.plot.png"),res=90,width=1000,height=1000)
sc <- autoEstCont(sc, doPlot=TRUE)
dev.off()
print("Infer corrected table of counts and rount to integer")
out <- adjustCounts(sc, roundToInt = TRUE)

print("save corrected count matrix")
scipy_sparse$save_npz(out_npz, out)

print("Done!")