#!/usr/bin/env Rscript

library("optparse")
library(HTSanalyzeR2)
library(org.Mm.eg.db)
library(KEGGREST)
library(igraph)
library(ggplotify)
library(GO.db)
library(limma)

# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("TxDb.Mmusculus.UCSC.mm10.knownGene")

library(TxDb.Mmusculus.UCSC.mm10.knownGene)

option_list = list(
  make_option(c("-gm", "--gene_matrix"), type="character", default=NULL, help="diann gene matrix output file name", metavar="character"),
  make_option(c("-fc", "--fc_col"), type="character", default=NULL, help="Sepcifify a column header with logfold to perform gsea", metavar="character"),
  make_option(c("-o", "--out"), type="character", default=NULL, help="output dir", metavar="character"),
  make_option(c("-sp", "--specie"), type="character", default=NULL, help="Species: Hs, Mm", metavar="character"),
  make_option(c("-gs", "--geneset"), type="character", default=NULL, help="Gene set: kegg, go, MSigDB", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$gene_matrix)){
  print_help(opt_parser)
  stop("script requires a diann gene matrix file", call.=FALSE)
}

f <- opt$gene_matrix

# test
f <- "/storage/Documents/service/externe/ilan/20240702_mouse_ms_organoid/gsean_outliersfiltered/all_folds.tsv"
gs <- "kegg"
gs <- "go"
gs <- "MSigDB"
# o <- paste("/storage/Documents/service/externe/ilan/20240702_mouse_ms_organoid/gsea/",gs,sep='')
o <- "/storage/Documents/service/externe/ilan/20240702_mouse_ms_organoid/gsean_outliersfiltered"
sp <- "Mm"
fc_col <- "fc"

g_list= read.csv(
  f, 
  header = TRUE,
  sep="\t",
  row.names = 1,
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
# ListGSC <- list(PW_KEGG=PW_KEGG)
# ListGSC <- list(GO=GO, PW_KEGG=PW_KEGG, MSig_C2=MSig_C2)

gsca <- GSCA(
  listOfGeneSetCollections=ListGSC, 
  geneList=phenotype
)
            
## preprocess
gsca1 <- preprocess(
  gsca, 
  species=sp, 
  initialIDs="SYMBOL",
  keepMultipleMappings=TRUE, 
  duplicateRemoverMethod="max",
  orderAbsValue=FALSE
)

## analysis
if (requireNamespace("doParallel", quietly=TRUE)) {
    doParallel::registerDoParallel(cores=4)
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
}



## draw GSEA plot for a specific gene set
out <- paste(o,gs,sep='/')
report(gsca3,reportDir = out, gseaPlot = FALSE)
t <- list.dirs(path = o, recursive = FALSE)
t <- grep(gs, t, value = TRUE)
plotGSEA(gsca3, gscs=gscs, ntop=10, filepath=t)



