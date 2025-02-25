
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
#
# BiocManager::install("edgeR")
# BiocManager::install("MAST")


suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("edgeR"))
# suppressPackageStartupMessages(library("MAST"))

fit_model <- function(c_data,g_data){
  # create an edgeR object with counts and grouping factor
  y <- DGEList(c_data, group = g_data$sample)
  # filter out genes with low counts
  print("Dimensions before subsetting:")
  print(dim(y))
  print("")
  keep <- filterByExpr(y)
  y <- y[keep, , keep.lib.sizes=FALSE]
  print("Dimensions after subsetting:")
  print(dim(y))
  print("")
  # normalize
  y <- calcNormFactors(y)
  # create a vector that is a concatenation of condition and cell type that we will later use with contrasts
  group <- paste0(g_data$sample, ".", g_data$cell_type)
  # replicate <- colData(adata_)$replicate
  # create a design matrix: here we have multiple donors so also consider that in the design matrix
  # design <- model.matrix(~ 0 + group + replicate)
  design <- model.matrix(~ 0 + group)
  # estimate dispersion
  y <- estimateDisp(y, design = design)
  # fit the model
  fit <- glmQLFit(y, design)
  return(list("fit"=fit, "design"=design, "y"=y))
}

option_list <- list(
  make_option(c("-o", "--outpath"), type="character", default=NULL, help="output file path", metavar="character"),
  make_option(c("-m", "--count_tsv"), type="character", default=NULL, help="Count tsv matrix", metavar="character"),
  make_option(c("-g", "--groups"), type="character", default=NULL, help="Groups tsv", metavar="character"),
  make_option(c("-r", "--res"), type="character", default=NULL, help="Leiden resolution", metavar="character"),
  make_option(c("-t", "--treated_regex"), type="character", default=NULL, help="treated group selection regex", metavar="character"),
  make_option(c("-c", "--ctrl_regex"), type="character", default=NULL, help="control group selection regex", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

out <- opt$outpath
counts <- opt$count_tsv
groups <- opt$groups
res <- opt$res
exp_regex <- opt$treated_regex
ctrl_regex <- opt$ctrl_regex

# counts <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/06-DifferentialGeneExpression/edger_counts.tsv"
# groups <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/06-DifferentialGeneExpression/edger_groups.tsv"
# out <- "/storage/Documents/service/externe/ilan/20241209_scRNAseq_FL_SA/out/06-DifferentialGeneExpression"
# exp_regex <- "FL"
# ctrl_regex <- "V"


print(paste0("out_path=",out))
print(paste0("in_counts=",counts))
print(paste0("in_groups=",groups))

print("loading counts")
c_data <- read.csv(
  counts,
  header = TRUE,
  row.names = 1,
  sep="\t",
  na.strings = "",
  stringsAsFactors=FALSE
)

print("loading groups list")
g_data <- read.csv(
  groups,
  header = TRUE,
  row.names = 1,
  sep="\t",
  na.strings = "",
  stringsAsFactors=FALSE
)

print("Fitting model")
model <- fit_model(c_data,g_data)

fit <- model$fit
y <- model$y
d <- model$design

pdf(paste0(out, "/Edger.plots.",res,".pdf"))
plotMDS(y, col=ifelse(y$samples$group == "FL_T", "red", "blue"))
plotBCV(y)
dev.off()

design_names <- colnames(y$design)
for (cell_type in unique(g_data$cell_type)) {
  print(cell_type)
  print(paste0("### Cell type: ", cell_type, " ###"))
  exp_names <- design_names[grepl(paste0("^group",exp_regex,"_.*",cell_type,"$"), colnames(y$design), perl = TRUE)] #No need for the * as you want anything that contains "FL_T"
  ctrl_names <- design_names[grepl(paste0("^group",ctrl_regex,"_.*",cell_type,"$"), colnames(y$design), perl = TRUE)] #No need for the * as you want anything that contains "V_T"
  print(paste0("Experiment group names:",exp_names))
  print(paste0("Control group names:",ctrl_names))
  
  for (e in exp_names) {
    for (c in ctrl_names) {
      # create contrast for this cell type
      contrast_name <- paste0(e, "-", c)
      print(paste0("Running contrast:",contrast_name))
      myContrast <- makeContrasts(contrast_name, levels = y$design)
      # perform QLF test
      qlf <- glmQLFTest(fit, contrast=myContrast)
      # get all of the DE genes and calculate Benjamini-Hochberg adjusted FDR
      tt <- topTags(qlf, n = Inf)
      
      # save in the list with the results for all the cell types
      exp_parts <- strsplit(e, "\\.")[[1]]
      exp_lbl <- gsub("^group(.*)$", "\\1", exp_parts, , perl = TRUE)
      ctrl_parts <- strsplit(c, "\\.")[[1]]
      ctrl_lbl <- gsub("^group(.*)$", "\\1", ctrl_parts, , perl = TRUE)
      
      write.table(
        tt$table,
        paste0(out,'/dge_', exp_lbl[1], "--", ctrl_lbl[1], "--", cell_type, "--", res , '.tsv'),
        sep = '\t', col.names=NA
      )
    }
  }
}

print("Done!")