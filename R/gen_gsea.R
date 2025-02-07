#!/usr/bin/env Rscript

suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library(HTSanalyzeR2))
suppressPackageStartupMessages(library(org.Mm.eg.db))
suppressPackageStartupMessages(library(KEGGREST))
suppressPackageStartupMessages(library(igraph))
suppressPackageStartupMessages(library(ggplotify))
suppressPackageStartupMessages(library(GO.db))
suppressPackageStartupMessages(library(limma))

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

# test
# f <- "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take3/EPostvsEPre.stats.tsv"
# gs <- "kegg"
# # gs <- "go"
# # gs <- "MSigDB"
# # o <- paste("/storage/Documents/service/externe/ilan/20240702_mouse_ms_organoid/gsea/",gs,sep='')
# o <- "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take3/gsea"
# sp <- "Hs"
# fc_col <- "log2fc"
# lbl <- "EPostvsEPre"

if (fc_col == ''){
    fc_col <- paste(lbl,"_diff",sep = "")
}

if (sp == 'Mm'){
  library(TxDb.Mmusculus.UCSC.mm10.knownGene)
} else if (sp == 'Hs'){
  library(TxDb.Hsapiens.UCSC.hg38.knownGene)
} else {
  message (paste0("Unrecongnised sp provided. Supported species are: Hs and Mm"))
  quit(1)
}

out <- paste(
    o,
    lbl,
    gs,
    sep='/'
)
if (!dir.exists(out)){
  dir.create(out, showWarnings = TRUE, recursive = TRUE)
}

g_list = read.csv(
  f, 
  header = TRUE,
  sep="\t",
  stringsAsFactors=FALSE
)

# Function to find the name of specific column labels
find_column_name_and_index <- function(df, labels) {
  for (label in labels) {
    index <- which(colnames(df) == label)
    if (length(index) > 0) {
      return(list(name = label, index = index))
    }
  }
  return(NULL)
}

# Define the labels to search for
labels_to_find <- c("UniProt", "Protein.Group")
# Get the name of the column
column_info <- find_column_name_and_index(g_list, labels_to_find)

# Check the column information
print(paste("column_name is:", column_info$name))
print(paste("column_index is:", column_info$index))

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
names(phenotype) <- as.vector(g_list[[column_info$name]])

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




