library(pathfindR)
library(dplyr)

gen_reactome <- function(f, o, fc_col, pval_col, gs_source, sp) {
  g_list= read.csv(
    f, 
    header = TRUE,
    sep="\t",
    na.strings = "", 
    stringsAsFactors=FALSE
  )
  
  g_list <- g_list %>% select(gene, fc_col, pval_col) %>% filter_at(vars(pval_col), all_vars(!is.na(.)))
  
  gsets_list <- get_gene_sets_list(
    source = gs_source,
    species = sp
  )
  
  output_df <- run_pathfindR(
    g_list,
    output_dir = o,
    custom_genes = gsets_list$gene_sets,
    custom_descriptions = gsets_list$descriptions,
    min_gset_size = 10,
    max_gset_size = 300,
    n_processes = 4,
    #pin_name_path = "STRING"
    pin_name_path = "Biogrid"
  )
  
  output_df_clustered <- cluster_enriched_terms(output_df, plot_dend = FALSE, plot_clusters_graph = FALSE)
  out=paste(o,"/enrichment_chart.tsv",sep='')
  write.table(output_df_clustered, file = out, sep = "\t")
  
  # plotting only selected clusters for better visualization
  selected_clusters <- subset(output_df_clustered[output_df_clustered$Status == "Representative", ], Cluster %in% 1:10)
  
  # output png enrichment chart
  out=paste(o,"/enrichment_chart.png",sep='')
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
  out=paste(o,"/enrichment_chart.svg",sep='')
  svg(
    filename = out,
    width = 8,
    height = 4
  )
  plot(enrichment_chart(selected_clusters, plot_by_cluster = TRUE))
  dev.off()
}

### analysis

# common param
f <- "/storage/Documents/service/externe/ilan/20230606_Mouse_DIA_MS/reactome/gene_list.tsv"
fc_col <- "NLRC5d_fc"
pval_col <- "ttest_NLRC5dvsWT"
sp <- "Mus musculus"

o <- "/storage/Documents/service/externe/ilan/20230606_Mouse_DIA_MS/reactome"
# Available gene sets are 'KEGG', 'Reactome', 'BioCarta', 'GO-All', 'GO-BP', 'GO-CC', 'GO-MF', 'cell_markers', 'mmu_KEGG'
for (gs_source in c('KEGG', 'Reactome', 'BioCarta', 'GO-All', 'GO-BP', 'GO-CC', 'GO-MF', 'cell_markers', 'mmu_KEGG')) {
  out <- paste(o,gs_source,sep='/')
  print(paste("##### Running",gs_source, sep = " "))
  gen_reactome(f, out, fc_col, pval_col, gs_source, sp)
}



