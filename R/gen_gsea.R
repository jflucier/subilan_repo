#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(HTSanalyzeR2))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(limma))

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")

library(TxDb.Mmusculus.UCSC.mm10.knownGene)

option_list = list(
  make_option(c("-m", "--gene_matrix"), type="character", default=NULL, help="diann gene matrix output file name", metavar="character"),
  make_option(c("-c", "--comp_label"), type="character", default="", help="Specify comparison label. Used to set pval_col and fc_col if not provided (i.e. <comp_label>_diff)", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="output dir", metavar="character"),
  make_option(c("-s", "--specie"), type="character", default=NULL, help="Species: Hs, Mm", metavar="character"),
  make_option(c("-f", "--fc_col"), type="character", default="", help="Sepcifify a column header with log2fc", metavar="character"),
  make_option(c("-g", "--geneset"), type="character", default=NULL, help="Gene set: kegg, go, MSigDB", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$gene_matrix)){
  print_help(opt_parser)
  stop("script requires a diann gene matrix file", call.=FALSE)
}

f <- opt$gene_matrix
lbl <- opt$comp_label
o <- opt$out
sp <- opt$specie
gs <- opt$geneset
fc_col <- opt$fc_col

if (fc_col == ''){
    fc_col <- paste(lbl,"_diff",sep = "")
}

out <- paste(
    o,
    comp_label,
    sep='/'
)
if (!dir.exists(out)){
  dir.create(out, showWarnings = TRUE, recursive = TRUE)
}

# test
# f <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/report.pg_matrix.proteoptypic.dge.tsv"
# gs <- "kegg"
# # gs <- "go"
# # gs <- "MSigDB"
# # o <- paste("/storage/Documents/service/externe/ilan/20240702_mouse_ms_organoid/gsea/",gs,sep='')
# o <- "/storage/Documents/service/externe/sheela/20240729_mouse_ms_lysM/results_rmoutliers/gsea"
# sp <- "Mm"
# fc_col <- "X15KO_vs_WT_diff"

g_list= read.csv(
  f, 
  header = TRUE,
  sep="\t",
  row.names = 2,
  stringsAsFactors=FALSE
)

## prepare input for analysis
if (gs == "kegg") {
  message ("loading kegg geneset")
  PW_KEGG <- KeggGeneSets(species=sp)
  ListGSC <- list(PW_KEGG=PW_KEGG)
  gscs <- c("PW_KEGG")
} else if (gs == "go") {
  message ("loading go geneset")
  GO <- GOGeneSets(species=sp, ontologies = c("BP","MF"))
  ListGSC <- list(GO=GO)
  gscs <- c("GO")
} else {
  message ("loading MSigDB geneset")
  MSig_C2 <- MSigDBGeneSets(collection = "C2", species = sp)
  ListGSC <- list(MSig_C2=MSig_C2)
  gscs=c("MSig_C2")
}

phenotype <- as.vector(g_list[[fc_col]])
names(phenotype) <- rownames(g_list)

gsca <- GSCA(
  listOfGeneSetCollections=ListGSC, 
  geneList=phenotype
)
            
## preprocess
gsca1 <- preprocess(
  gsca, 
  species=sp, 
  initialIDs="UNIPROT",
  keepMultipleMappings=TRUE, 
  duplicateRemoverMethod="max",
  orderAbsValue=FALSE
)

## analysis
if (requireNamespace("doParallel", quietly=TRUE)) {
    doParallel::registerDoParallel(cores=8)
}  ## support parallel calculation using multiple cores
gsca2 <- analyze(
  gsca1, 
  para=list(
      pValueCutoff=0.05, 
      pAdjustMethod="BH",
      nPermutations=1000, 
      minGeneSetSize=5,
      exponent=1),
  doGSOA = FALSE
)

## append gene sets terms
if (gs == "kegg") {
  message ("append kegg gene sets terms")
  gsca3 <- appendGSTerms(
    gsca2, 
    keggGSCs=gscs,
    species=sp
  )
  
  topGS <- getTopGeneSets(
    gsca3, 
    resultName="GSEA.results",
    gscs=gscs,
    allSig=TRUE
  )
  
  d <- gsca3@result[["GSEA.results"]][["PW_KEGG"]]
} else if (gs == "go") {
  message ("append go gene sets terms")
  gsca3 <- appendGSTerms(
    gsca2, 
    goGSCs=gscs,
    species=sp
  )
  
  topGS <- getTopGeneSets(
    gsca3, 
    resultName="GSEA.results",
    gscs=gscs,
    allSig=TRUE
  )
  d <- gsca3@result[["GSEA.results"]][["GO"]]
} else {
  message ("append MSig gene sets terms")
  gsca3 <- appendGSTerms(
    gsca2, 
    msigdbGSCs = gscs,
    species=sp
  )
  
  topGS <- getTopGeneSets(
    gsca3, 
    resultName="GSEA.results",
    gscs=gscs,
    allSig=TRUE
  )
  d <- gsca3@result[["GSEA.results"]][["MSig_C2"]]
}



## draw GSEA plot for a specific gene set
print("Outputting GSEA")
if (length(topGS) > 0) {
  plotGSEA(gsca3, gscs=gscs, filepath=out, allSig=TRUE)
}

print("Outputting geneset TSV")
write.table(
  data.frame("Gene.Set"=rownames(d),d),
  file = paste0(out,"/genesets.tsv"),
  sep = "\t",
  row.names = FALSE
)

print("done!")

# out <- paste(o,n,sep='/')
# report(gsca3,reportDir = out, gseaPlot = TRUE)
# summarize(gsca3)




