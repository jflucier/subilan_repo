# install.packages("SomaDataIO")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE)) {
#   install.packages("BiocManager")
# }
# BiocManager::install("Biobase")
# 
# install.packages("tidyverse")
suppressPackageStartupMessages(library("optparse"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("SomaDataIO"))

###### Output adat data to tsv #####
option_list <- list(
  make_option(c("-o", "--out"), type="character", default=NULL, help="output file path", metavar="character"),
  make_option(c("-a", "--adat"), type="character", default=NULL, help="Input files filtered data matrix npz file", metavar="character"),
  make_option(c("-s", "--samples"), type="character", default=NULL, help="Input files raw data matrix npz file", metavar="character"),
  make_option(c("-c", "--comp_name"), type="character", default=NULL, help="Group column label", metavar="character")
);

opt_parser <- OptionParser(option_list=option_list);
opt <- parse_args(opt_parser);

outpath <- opt$out
adat_f <- opt$adat
annot_f <- opt$samples
comp_name <- opt$comp_name

# adat_f <- "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take3/SS-2453319_v4.1_EDTAPlasma.hybNorm.medNormInt.plateScale.calibrate.anmlQC.qcCheck.anmlSMP.20240528.adat"
# annot_f <- "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take3/plate_sample_annot.tsv"
# outpath <- "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take3"

my_adat <- read_adat(adat_f)
annot_tbl <- read.csv(
  annot_f, 
  header = TRUE,
  sep="\t",
  row.names = 1,
  stringsAsFactors=FALSE
)

my_adat[ , 'SubjectID_new'] <- NA
my_adat$SubjectID_new <- annot_tbl[match(row.names(my_adat), row.names(annot_tbl)),"SubjectID_new"]

my_adat[ , 'group'] = NA
my_adat$group <- annot_tbl[match(row.names(my_adat), row.names(annot_tbl)),"group"]

my_adat[ , 'sex'] = NA
my_adat$sex <- annot_tbl[match(row.names(my_adat), row.names(annot_tbl)),"sex"]

my_adat[ , 'exercice'] = NA
my_adat$exercice <- annot_tbl[match(row.names(my_adat), row.names(annot_tbl)),"exercice"]

my_adat[ , 'CompE'] = NA
my_adat$CompE <- annot_tbl[match(row.names(my_adat), row.names(annot_tbl)),"CompE"]

my_adat[ , 'CompNE'] = NA
my_adat$CompNE <- annot_tbl[match(row.names(my_adat), row.names(annot_tbl)),"CompNE"]

cleanData <- my_adat |>
  filter(SampleType == "Sample") |>
  drop_na((!!as.symbol(comp_name)))

raw <- as.data.frame(cleanData)
# output raw data to tsv
out=paste(outpath,"/adat.raw.t.tsv",sep='')
write.table(
  raw,
  file = out,
  sep = "\t",
  row.names = TRUE,
  col.names = TRUE
)

# output target genes
target_map_all <-  getAnalyteInfo(cleanData) |>
  select(AptName, SeqId, Target = TargetFullName, EntrezGeneSymbol, UniProt)
tm <- as.data.frame(target_map_all)
out=paste(outpath,"/all_targets.tsv",sep='')
write.table(
  tm,
  file = out,
  sep = "\t",
  row.names = FALSE,
  col.names = TRUE
)

