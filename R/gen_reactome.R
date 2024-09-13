library(pathfindR)
library(dplyr)
library(tidyverse)
library("optparse")

gen_reactome <- function(f, o, gene_col, fc_col, pval_col, gs_source, sp_lbl, pin) {
  g_list= read.csv(
    f, 
    header = TRUE,
    sep="\t",
    na.strings = "NA", 
    stringsAsFactors=FALSE
  )
  g_list <- g_list %>% 
    select(all_of(gene_col), all_of(fc_col), all_of(pval_col)) %>%
    filter_at(vars(pval_col), all_vars(!is.na(.)))
  
  print(paste("Gene source is ", gs_source, sep= ""))
  if (gs_source == "KEGG") {
    if(sp_lbl == "human"){
      sp <- 'hsa'
    } else if (sp_lbl == "mouse") {
      sp <- 'mmu'
    } else{
      stop(paste("Unrecongnised specie provided. Possible values are: human or mouse. Value passed: ",sp_lbl,sep=''))
    }
    
    message(paste("Gene set list set for KEGG with organism code ", sp, sep= ""))
    gsets_list <- get_gene_sets_list(
      source = gs_source,
      org_code = sp
    )
  } else if (gs_source == "Reactome") {
    message(paste("Gene set list set for Reactome", sp, sep= ""))
    gsets_list <- get_gene_sets_list(
      source = gs_source
    )
  } else if (gs_source == "MSigDB") {
    if(sp_lbl != "human" && sp_lbl != "mouse"){
      stop(paste("Unrecongnised specie provided. Possible values are: human or mouse. Value passed: ",sp_lbl,sep=''))
    }
    message(paste("Gene set list set for MSigDB using collection C2 for specie ", sp, sep= ""))
    gsets_list <- get_gene_sets_list(
      source = gs_source,
      species = sp_lbl,
      collection = "C2"
    )
  }
  
  dir.create(o, showWarnings = TRUE, recursive = TRUE)
  d <- paste(o,pin,sep='/')
  # dir.create(d, showWarnings = TRUE, recursive = TRUE)
  output_df <- run_pathfindR(
    g_list,
    output_dir = d,
    custom_genes = gsets_list$gene_sets,
    custom_descriptions = gsets_list$descriptions,
    min_gset_size = 10,
    max_gset_size = 300,
    n_processes = 4,
    #pin_name_path = "STRING"
    pin_name_path = pin
  )
  
  output_df_clustered <- cluster_enriched_terms(output_df, plot_dend = FALSE, plot_clusters_graph = FALSE)
  out=paste(d,"enrichment_chart.tsv",sep='/')
  write.table(output_df_clustered, file = out, sep = "\t")
  
  # plotting only selected clusters for better visualization
  selected_clusters <- subset(output_df_clustered[output_df_clustered$Status == "Representative", ], Cluster %in% 1:10)
  
  # output png enrichment chart
  out=paste(d,"enrichment_chart.png",sep='/')
  png(
    filename = out,
    res = 250,
    width = 8,
    height = 4,
    units = "in"
  )
  plot(enrichment_chart(selected_clusters, plot_by_cluster = TRUE))
  dev.off()
  
  # output png enrichment chart
  out=paste(d,"enrichment_chart.svg",sep='/')
  svg(
    filename = out,
    width = 8,
    height = 4
  )
  plot(enrichment_chart(selected_clusters, plot_by_cluster = TRUE))
  dev.off()
}
  #for (pin in c('Biogrid', 'STRING', 'GeneMania', 'IntAct', 'KEGG', 'mmu_STRING')) {
  # for (pin in c('Biogrid', 'STRING', 'IntAct', 'KEGG', 'mmu_STRING')) {
  #   print(paste("pin=",pin, sep= " "))
  #   d <- paste(o,pin,sep='/')
  #   # dir.create(d, showWarnings = TRUE, recursive = TRUE)
  #   output_df <- run_pathfindR(
  #     g_list,
  #     output_dir = d,
  #     custom_genes = gsets_list$gene_sets,
  #     custom_descriptions = gsets_list$descriptions,
  #     min_gset_size = 10,
  #     max_gset_size = 300,
  #     n_processes = 4,
  #     #pin_name_path = "STRING"
  #     pin_name_path = pin
  #   )
  # 
  #   output_df_clustered <- cluster_enriched_terms(output_df, plot_dend = FALSE, plot_clusters_graph = FALSE)
  #   out=paste(d,"enrichment_chart.tsv",sep='/')
  #   write.table(output_df_clustered, file = out, sep = "\t")
  #   
  #   # plotting only selected clusters for better visualization
  #   selected_clusters <- subset(output_df_clustered[output_df_clustered$Status == "Representative", ], Cluster %in% 1:10)
  #   
  #   # output png enrichment chart
  #   out=paste(d,"enrichment_chart.png",sep='/')
  #   png(
  #     filename = out,
  #     res = 250,
  #     width = 8,
  #     height = 4,
  #     units = "in"
  #   )
  #   plot(enrichment_chart(selected_clusters, plot_by_cluster = TRUE))
  #   dev.off()
  #   
  #   # output png enrichment chart
  #   out=paste(d,"enrichment_chart.svg",sep='/')
  #   svg(
  #     filename = out,
  #     width = 8,
  #     height = 4
  #   )
  #   plot(enrichment_chart(selected_clusters, plot_by_cluster = TRUE))
  #   dev.off()
  # }
# }


# f <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/report.pg_matrix.proteoptypic.dge.modif.tsv"
# sp <- "mouse"
# sp_lbl <- sp
# out <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/reactome/LysM_minus_vs_WT"
# gene_col <- "Genes"
# fc_col <- "LysM_minus_vs_WT_diff"
# pval_col <- "LysM_minus_vs_WT_p.adj"

option_list = list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="frapipe proteotypic protein group matrix", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="output dir", metavar="character"),
  make_option(c("-s", "--specie"), type="character", default=NULL, help="Species: human, mouse", metavar="character"),
  make_option(c("-g", "--gene_col"), type="character", default=NULL, help="Gene header column", metavar="character"),
  make_option(c("-f", "--fc_col"), type="character", default=NULL, help="Sepcifify a column header with logfold", metavar="character"),
  make_option(c("-p", "--pval_col"), type="character", default=NULL, help="Sepcifify a column header with pvalues", metavar="character"),
  make_option(c("--gene_source"), type="character", default=NULL, help="Gene source: KEGG, Reactome, MSigDB", metavar="character"),
  make_option(c("--pin"), type="character", default=NULL, help="PIN: Biogrid, STRING, IntAct, KEGG, mmu_STRING", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$input)){
  print_help(opt_parser)
  stop("script requires a diann input matrix file", call.=FALSE)
}

f <- opt$input
out <- opt$out
gene_col <- opt$gene_col
fc_col <- opt$fc_col
pval_col <- opt$pval_col
gs_source <- opt$gene_source
sp <- opt$specie
pin <- opt$pin

print(paste("##### Running",gs_source, sep = " "))
if (!dir.exists(out)){
  dir.create(out, showWarnings = TRUE, recursive = TRUE)
}
o <- paste(out,gs_source,sep='/')
gen_reactome(f, o, gene_col, fc_col, pval_col, gs_source, sp, pin)

# f <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/report.pg_matrix.proteoptypic.dge.modif.tsv"
# sp <- "mouse"
# sp_lbl <- sp
# out <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/reactome/X15KO_vs_WT"
# gene_col <- "Genes"
# fc_col <- "X15KO_vs_WT_diff"
# pval_col <- "X15KO_vs_WT_p.adj"
# 
# dir.create(out, showWarnings = TRUE, recursive = TRUE)
# for (gs_source in c('KEGG', 'Reactome', 'MSigDB')) {
#   tmp <- out
#   o <- paste(tmp,gs_source,sep='/')
#   print(paste("##### Running",gs_source, sep = " "))
#   gen_reactome(f, o, gene_col, fc_col, pval_col, gs_source, sp)
# }
# 
# f <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/report.pg_matrix.proteoptypic.dge.modif.tsv"
# sp <- "mouse"
# sp_lbl <- sp
# out <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/reactome/LysM_minus_vs_WT"
# gene_col <- "Genes"
# fc_col <- "LysM_minus_vs_WT_diff"
# pval_col <- "LysM_minus_vs_WT_p.adj"
# 
# dir.create(out, showWarnings = TRUE, recursive = TRUE)
# for (gs_source in c('KEGG', 'Reactome', 'MSigDB')) {
#   tmp <- out
#   o <- paste(tmp,gs_source,sep='/')
#   print(paste("##### Running",gs_source, sep = " "))
#   gen_reactome(f, o, gene_col, fc_col, pval_col, gs_source, sp)
# }
# 
# f <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/report.pg_matrix.proteoptypic.dge.modif.tsv"
# sp <- "mouse"
# sp_lbl <- sp
# out <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/reactome/LysM_plus_vs_WT"
# gene_col <- "Genes"
# fc_col <- "LysM_plus_vs_WT_diff"
# pval_col <- "LysM_plus_vs_WT_p.adj"
# 
# dir.create(out, showWarnings = TRUE, recursive = TRUE)
# for (gs_source in c('KEGG', 'Reactome', 'MSigDB')) {
#   tmp <- out
#   o <- paste(tmp,gs_source,sep='/')
#   print(paste("##### Running",gs_source, sep = " "))
#   gen_reactome(f, o, gene_col, fc_col, pval_col, gs_source, sp)
# }
# 
# f <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/report.pg_matrix.proteoptypic.dge.modif.tsv"
# sp <- "mouse"
# sp_lbl <- sp
# out <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/reactome/X15KO_vs_LysM_minus"
# gene_col <- "Genes"
# fc_col <- "X15KO_vs_LysM_minus_diff"
# pval_col <- "X15KO_vs_LysM_minus_p.adj"
# 
# dir.create(out, showWarnings = TRUE, recursive = TRUE)
# for (gs_source in c('KEGG', 'Reactome', 'MSigDB')) {
#   tmp <- out
#   o <- paste(tmp,gs_source,sep='/')
#   print(paste("##### Running",gs_source, sep = " "))
#   gen_reactome(f, o, gene_col, fc_col, pval_col, gs_source, sp)
# }
# 
# f <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/report.pg_matrix.proteoptypic.dge.modif.tsv"
# sp <- "mouse"
# sp_lbl <- sp
# out <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/reactome/X15KO_vs_LysM_plus"
# gene_col <- "Genes"
# fc_col <- "X15KO_vs_LysM_plus_diff"
# pval_col <- "X15KO_vs_LysM_plus_p.adj"
# 
# dir.create(out, showWarnings = TRUE, recursive = TRUE)
# for (gs_source in c('KEGG', 'Reactome', 'MSigDB')) {
#   tmp <- out
#   o <- paste(tmp,gs_source,sep='/')
#   print(paste("##### Running",gs_source, sep = " "))
#   gen_reactome(f, o, gene_col, fc_col, pval_col, gs_source, sp)
# }


### analysis

# common param
# f <- "/storage/Documents/service/externe/ilan/20230606_Mouse_DIA_MS/reactome/gene_list.tsv"
# # specie can be human or mouse
# sp <- "mouse"
# sp_lbl <- sp

### NLRC5d
# o <- "/storage/Documents/service/externe/ilan/20230606_Mouse_DIA_MS/reactome/NLRC5d"
# fc_col <- "NLRC5d_fc"
# pval_col <- "ttest_NLRC5dvsWT"
# 
# gs_source <- 'KEGG'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
# gs_source <- 'Reactome'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
# gs_source <- 'MSigDB'
# oo <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, oo, fc_col, pval_col, gs_source, sp)

### RAG1d
# o <- "/storage/Documents/service/externe/ilan/20230606_Mouse_DIA_MS/reactome/RAG1d"
# fc_col <- "RAG1d_fc"
# pval_col <- "ttest_RAG1dvsWT"
# 
# gs_source <- 'KEGG'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
# gs_source <- 'Reactome'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
# gs_source <- 'MSigDB'
# oo <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, oo, fc_col, pval_col, gs_source, sp)


# ### NLRC5 vs RAG1d
# f <- "/storage/Documents/service/externe/ilan/20230606_Mouse_DIA_MS/reactome/gene_list_nlrc5vsrag1.tsv"
# # specie can be human or mouse
# sp <- "mouse"
# sp_lbl <- sp
# o <- "/storage/Documents/service/externe/ilan/20230606_Mouse_DIA_MS/reactome/NLRC5vsRAG1d"
# fc_col <- "fc_NLRC5dvsRAG1d"
# pval_col <- "ttest_NLRC5dvsRAG1d"
# 
# gs_source <- 'KEGG'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
# gs_source <- 'Reactome'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
# gs_source <- 'MSigDB'
# oo <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, oo, fc_col, pval_col, gs_source, sp)

# ### 20240702_mouse_ms_organoid
# f <- "/storage/Documents/service/externe/ilan/20240702_mouse_ms_organoid/reactome/gene_list.tsv"
# # specie can be human or mouse
# sp <- "mouse"
# sp_lbl <- sp
# o <- "/storage/Documents/service/externe/ilan/20240702_mouse_ms_organoid/reactome"
# fc_col <- "fc_KOvsWT"
# pval_col <- "ttest_KOvsWT"
# 
# gs_source <- 'KEGG'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
# gs_source <- 'Reactome'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
# gs_source <- 'MSigDB'
# oo <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, oo, fc_col, pval_col, gs_source, sp)
# 
# ### 20240729_mouse_ms_lysM
# f <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/report.pg_matrix.proteoptypic.dge.modif.tsv"
# # specie can be human or mouse
# sp <- "mouse"
# sp_lbl <- sp
# o <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/reactome"
# gene_col <- "Genes"
# fc_col <- "X15KO_vs_WT_diff"
# pval_col <- "X15KO_vs_WT_p.adj"
# 
# gs_source <- 'KEGG'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, gene_col, fc_col, pval_col, gs_source, sp)
# gs_source <- 'Reactome'
# out <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
# gs_source <- 'MSigDB'
# oo <- paste(o,gs_source,sep='/')
# print(paste("##### Running",gs_source, sep = " "))
# gen_reactome(f, oo, fc_col, pval_col, gs_source, sp)


