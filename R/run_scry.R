
# scry: Small-Count Analysis Methods for High-Dimensional Data
# https://www.bioconductor.org/packages/release/bioc/html/scry.html

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# BiocManager::install("scry")
# BiocManager::install("SummarizedExperiment")

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("scry"))
suppressPackageStartupMessages(library("SummarizedExperiment"))
suppressPackageStartupMessages(library("reticulate"))

sc <- import("scanpy")
np <- import("numpy")

option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default=NULL, help="output file path", metavar="character"),
  make_option(c("-i", "--in_h5ad"), type="character", default=NULL, help="Input h5ad files representing normalised adata object", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

out <- opt$outpath
in_h5ad <- opt$in_h5ad

# in_h5ad <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/adata_normalisation.h5ad"
# out <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/FeatureSelection"

print(paste0("out_path=",out))
print(paste0("in_h5ad=",in_h5ad))

print("Reading H5AD file")
adata <- sc$read_h5ad(in_h5ad)
print("Running devianceFeatureSelection")
binomial_deviance <- devianceFeatureSelection(adata$T$X)

# idx <- tail(np$argsort(binomial_deviance),4000)

idx <- tail(sort(binomial_deviance, index.return = TRUE)$ix,2000)

mask = np$zeros(adata$var_names$shape, dtype=np$bool_)
mask[idx] = TRUE

print(paste0("Outputting results in ",out))
write.table(
  mask,
  file = paste0(out,"/mask.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

write.table(
  binomial_deviance,
  file = paste0(out,"/binomial_deviance.tsv"),
  sep = "\t",
  row.names = FALSE,
  col.names = FALSE
)

print("Done!")