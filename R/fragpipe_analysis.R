
# FragPipeAnalystR to install and test
# https://www.nesvilab.org/FragPipeAnalystR/global_DIA_prot_tutorial.html

# install.packages("renv")
# 
# renv::install("bioc::ComplexHeatmap", prompt=FALSE)
# renv::install("bioc::limma", prompt=FALSE)
# renv::install("bioc::MSnbase", prompt=FALSE)
# renv::install("bioc::SummarizedExperiment", prompt=FALSE)
# renv::install("bioc::cmapR", prompt=FALSE)
# renv::install("bioc::ConsensusClusterPlus", prompt=FALSE)
# renv::install("Nesvilab/FragPipeAnalystR", prompt=FALSE)
# 
# # optional
# renv::install("nicolerg/ssGSEA2", prompt=FALSE)

library(FragPipeAnalystR)
library("optparse")

### first filter report.pg_matrix.tsv to select proteotypic proteins groups ###
# perl -ne '
# chomp($_);
# my @t = split("\t",$_);
# my @prot_ident = split(";",$t[0]);
# if(scalar(@prot_ident) == 1){
#   print $_ . "\n";
# }
# ' report.pg_matrix.tsv > report.pg_matrix.proteoptypic.tsv

# bp <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM"
# out_bp <- paste0(bp,"/fragpipe")
# in_diann <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results/report.pg_matrix.proteoptypic.tsv"
# in_design <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/experiment_annotation.tsv"

option_list = list(
  make_option(c("-o", "--out"), type="character", default=NULL, help="output dir", metavar="character"),
  make_option(c("-m", "--matrix"), type="character", default=NULL, help="diann pg matrix filtered for proteotypic proteins", metavar="character"),
  make_option(c("-d", "--design"), type="character", default=NULL, help="Fragpipe TSV design file. Header is: file\tsample\tsample_name\tcondition\treplicate", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

out_bp <- opt$out
in_diann <- opt$matrix
in_design <- opt$design

if (!dir.exists(out_bp)){
  dir.create(out_bp, showWarnings = TRUE)
}

ccrcc <- make_se_from_files(
  in_diann,
  in_design,
  type = "DIA",
  level = "protein"
)

pdf(
  paste0(out_bp,"/report.pg_matrix.proteoptypic.qc.pdf")
)
plot_pca(
  ccrcc
)

plot_correlation_heatmap(ccrcc)
plot_missval_heatmap(ccrcc)
plot_feature_numbers(ccrcc)
dev.off()

# manual_impute imputes missing values in a proteomics dataset by random draws from a manually defined distribution.
pdf(
  paste0(out_bp,"/report.pg_matrix.proteoptypic.qc_imputed.pdf")
)
imputed <- manual_impute(ccrcc)
plot_pca(
  imputed
)

plot_correlation_heatmap(imputed)
plot_feature_numbers(imputed)
dev.off()

pdf(
  paste0(out_bp,"/report.pg_matrix.proteoptypic.dge.pdf")
)
de_result <- test_limma(ccrcc, type = "all")
# add_rejections marks significant proteins based on defined cutoffs. Deafult is <= 0.05 adjusted P value and log2 fold change >= 1
de_result_updated <- add_rejections(de_result)

group_list <- unique(de_result@colData@listData[["condition"]])
for (value1 in group_list) {
  #print(value1)
  for (value2 in group_list) {
    v_diff <- paste0(value1,"_vs_",value2,"_diff")
    if(!is.null(de_result_updated@elementMetadata@listData[[v_diff]])){
      lbl <- paste0(value1,"_vs_",value2)
      print(paste0("plotting volcano for: ", lbl))
      plot_volcano(de_result_updated, lbl, name_col = "Genes")
    }
  }
}
# plot_volcano(de_result_updated, "X15KO_vs_WT", name_col = "Genes")
# plot_volcano(de_result_updated, "LysM_minus_vs_WT", name_col = "Genes")
# plot_volcano(de_result_updated, "LysM_plus_vs_WT", name_col = "Genes")
# plot_volcano(de_result_updated, "LysM_minus_vs_LysM_plus", name_col = "Genes")
# plot_volcano(de_result_updated, "X15KO_vs_LysM_minus", name_col = "Genes")
# plot_volcano(de_result_updated, "X15KO_vs_LysM_plus", name_col = "Genes")
dev.off()

x <- de_result_updated@elementMetadata@listData
y <- data.frame("rowid"=x[["Protein.Group"]],x)
write.table(
  data.frame("rowid"=rownames(y),x),
  file = paste0(out_bp,"/report.pg_matrix.proteoptypic.dge.tsv"),
  sep = "\t",
  row.names = FALSE
)

pdf(
  paste0(out_bp,"/report.pg_matrix.proteoptypic.dge_imputed.pdf")
)
de_result <- test_limma(imputed, type = "all")
# add_rejections marks significant proteins based on defined cutoffs. Deafult is <= 0.05 adjusted P value and log2 fold change >= 1
de_result_updated <- add_rejections(de_result)
group_list <- unique(de_result@colData@listData[["condition"]])
for (value1 in group_list) {
  #print(value1)
  for (value2 in group_list) {
    v_diff <- paste0(value1,"_vs_",value2,"_diff")
    if(!is.null(de_result_updated@elementMetadata@listData[[v_diff]])){
      lbl <- paste0(value1,"_vs_",value2)
      print(paste0("plotting volcano for: ", lbl))
      plot_volcano(de_result_updated, lbl, name_col = "Genes")
    }
  }
}
# 
# plot_volcano(de_result_updated, "X15KO_vs_WT", name_col = "Genes")
# plot_volcano(de_result_updated, "LysM_minus_vs_WT", name_col = "Genes")
# plot_volcano(de_result_updated, "LysM_plus_vs_WT", name_col = "Genes")
# plot_volcano(de_result_updated, "LysM_minus_vs_LysM_plus", name_col = "Genes")
# plot_volcano(de_result_updated, "X15KO_vs_LysM_minus", name_col = "Genes")
# plot_volcano(de_result_updated, "X15KO_vs_LysM_plus", name_col = "Genes")
dev.off()

x <- de_result_updated@elementMetadata@listData
y <- data.frame("rowid"=x[["Protein.Group"]],x)
write.table(
  data.frame("rowid"=rownames(y),x),
  file = paste0(out_bp,"/report.pg_matrix.proteoptypic.dge_imputed.tsv"),
  sep = "\t",
  row.names = FALSE
)
