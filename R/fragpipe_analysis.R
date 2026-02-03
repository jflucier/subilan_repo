
# FragPipeAnalystR to install and test
# https://www.nesvilab.org/FragPipeAnalystR/global_DIA_prot_tutorial.html

# install.packages("renv")
# # 
# renv::install("bioc::ComplexHeatmap", prompt=FALSE)
# renv::install("bioc::limma", prompt=FALSE)
# renv::install("bioc::MSnbase", prompt=FALSE)
# renv::install("bioc::SummarizedExperiment", prompt=FALSE)
# renv::install("bioc::cmapR", prompt=FALSE)
# renv::install("bioc::ConsensusClusterPlus", prompt=FALSE)
# renv::install("Nesvilab/FragPipeAnalystR", prompt=FALSE)
# renv::install("bioc::SEtools", prompt=FALSE)
# 
# # optional
# renv::install("nicolerg/ssGSEA2", prompt=FALSE)
# install.packages("pals") 
# install.packages("Polychrome")

suppressPackageStartupMessages(library(FragPipeAnalystR))
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(SEtools))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(Polychrome))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(pals))

# bp <- "/storage/Documents/service/externe/sheela/20240829_HSCs_liver_mouse"
# out_bp <- paste0(bp,"/fragpipe")
# in_diann <- "/storage/Documents/service/externe/sheela/20240829_HSCs_liver_mouse/results/report.pg_matrix.proteoptypic.tsv"
# in_design <- "/storage/Documents/service/externe/sheela/20240829_HSCs_liver_mouse/experiment_annotation.outliers.tsv"

option_list = list(
  make_option(c("-f", "--format"), type="character", default="pdf", help="Output format: pdf or svg [default %default]", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="output dir", metavar="character"),
  make_option(c("-m", "--matrix"), type="character", default=NULL, help="diann pg matrix filtered for proteotypic proteins", metavar="character"),
  make_option(c("-d", "--design"), type="character", default=NULL, help="Fragpipe TSV design file. Header is: file\tsample\tsample_name\tcondition\treplicate", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

out_fmt <- opt$format
out_bp <- opt$out
in_diann <- opt$matrix
in_design <- opt$design

# out_bp <- "/storage/Documents/service/externe/sheela/20250107_mice_colon_take2/fragpipev1.0.1"
# in_diann <- "/storage/Documents/service/externe/sheela/20250107_mice_colon_take2/report.pg_matrix.tsv"
# in_design <- "/storage/Documents/service/externe/sheela/20250107_mice_colon_take2/experiment_annotation.tsv"

# Helper to open the correct device
open_device <- function(filename_base, format) {
  if (tolower(format) == "svg") {
    # %03d creates separate files for each plot (001, 002, etc.)
    svg(paste0(filename_base, "%03d.svg"), width = 10, height = 8)
  } else {
    pdf(paste0(filename_base, ".pdf"), width = 10, height = 8)
  }
}

output_dge <- function(condition, data, out) {
  group_list <- unique(condition)
  ctrl <- grep(
    "Ctrl", 
    group_list, 
    ignore.case=TRUE,
    value=TRUE, 
    fixed=FALSE
  )
  exp=grep(
    "Ctrl", 
    group_list, 
    ignore.case=TRUE,
    value=TRUE, 
    fixed=FALSE,
    invert=TRUE
  )
  
  i<-0
  for (c in group_list) {
    de_result_tmp <- test_limma(data, type = "control", control=c)
    if(i == 0){
      de_result <- de_result_tmp
    } else {
      de_result <- mergeSEs( list(se1=de_result, se2=de_result_tmp) )
    }
    i<-i+1
  }
  
  # add_rejections marks significant proteins based on defined cutoffs. Deafult is <= 0.05 adjusted P value and log2 fold change >= 1
  de_result_updated <- add_rejections(de_result)
  
  open_device(out, "pdf")
  for (value1 in group_list) {
    #print(value1)
    for (value2 in group_list) {
      v_diff <- paste0(value1,"_vs_",value2,"_diff")
      if(!is.null(de_result_updated@elementMetadata@listData[[v_diff]])){
        lbl <- paste0(value1,"_vs_",value2)
        print(paste0("plotting volcano for: ", lbl))
        p <- plot_volcano(de_result_updated, lbl, name_col = "Genes")
        print(p)
      }
    }
  }
  dev.off()
  
  x <- de_result_updated@elementMetadata@listData
  y <- data.frame("rowid"=x[["Protein.Group"]],x)
  write.table(
    data.frame("rowid"=rownames(y),x),
    file = paste0(out, ".tsv"),
    sep = "\t",
    row.names = FALSE
  )
  
}

if (!dir.exists(out_bp)){
  dir.create(out_bp, showWarnings = TRUE)
}

ccrcc <- make_se_from_files(
  in_diann,
  in_design,
  type = "DIA",
  level = "protein"
)

ccrcc_imputed <- manual_impute(ccrcc)
assay_data <- assay(ccrcc_imputed)

if ("batch" %in% colnames(colData(ccrcc_imputed))) {
  batch_info <- colData(ccrcc_imputed)$batch
  n_batches <- length(unique(batch_info))
  if(n_batches > 1){
    design_mat <- model.matrix(~condition, data=as.data.frame(colData(ccrcc_imputed)))
    assay(ccrcc_imputed) <- limma::removeBatchEffect(
      assay_data, 
      batch = batch_info, 
      design = design_mat
    )
  }
} 

condition_info <- colData(ccrcc_imputed)$condition
n_conditions <- length(unique(condition_info))

# 1. Define your custom palette (e.g., using the Polychrome built-in)
my_colors <- as.vector(palette.colors(n = n_conditions, palette = "Polychrome 36"))

# 2. Overwrite the default discrete scale function globally
# This works even if the wrapper function doesn't allow a palette argument
scale_colour_discrete <- function(...) scale_color_manual(values = my_colors)
scale_fill_discrete <- function(...) scale_fill_manual(values = my_colors)

open_device(paste0(out_bp,"/report.pg_matrix.proteoptypic.batchcorr.qc"), out_fmt)
plot_pca(ccrcc_imputed)
plot_correlation_heatmap(ccrcc_imputed)
plot_feature_numbers(ccrcc_imputed)
dev.off()

open_device(paste0(out_bp,"/report.pg_matrix.proteoptypic.qc"), out_fmt)
plot_pca(ccrcc)
plot_correlation_heatmap(ccrcc)
plot_missval_heatmap(ccrcc)
plot_feature_numbers(ccrcc)
dev.off()

rm(scale_colour_discrete, scale_fill_discrete)

outpath <- paste0(out_bp,"/report.pg_matrix.proteoptypic.dge")
output_dge(condition_info, ccrcc, outpath)

outpath <- paste0(out_bp,"/report.pg_matrix.proteoptypic.batchcorr.dge")
output_dge(condition_info, ccrcc_imputed, outpath)

