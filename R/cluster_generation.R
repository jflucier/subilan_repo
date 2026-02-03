
library(ggbiplot)
library(devtools)
library("RColorBrewer")
library(corrplot)
library(dplyr)

# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("ComplexHeatmap")

library("ComplexHeatmap")

# install.packages('gplots')
library("gplots")

all_files <- c(
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/E.adat.raw.diff.heatmap.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/E.adat.raw.heatmap.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/E.adat.raw.heatmap.filtered.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/E.adat.raw.diff.heatmap.filtered.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/NE.adat.raw.heatmap.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/NE.adat.raw.diff.heatmap.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/NE.adat.raw.heatmap.filtered.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/ALL.adat.raw.heatmap.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/ALL.adat.raw.heatmap.diff.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/ALL.adat.raw.heatmap.diff.filtered.tsv",
  "/storage/Documents/service/externe/sheela/20240606_LC_exercice_somalogic/results.take4/V2/ALL.adat.raw.heatmap.filtered.tsv"
)

# f <- all_files[7]

for (f in all_files) {
  message("running:",f)
  # 1. Read the file without row.names first to keep alignment raw
  all_data <- read.delim(f, header = FALSE, sep = "\t", skip = 1, check.names = FALSE)
  
  # 2. Extract the IDs from the first column and set them as row names
  rownames(all_data) <- all_data[, 1]
  
  # 3. Remove that first column so only numeric data remains
  all_data <- all_data[, -1]
  
  # 4. Manually re-assign the protein names from the file's first line
  # We read the first line and remove the first (empty) element
  header_row <- scan(f, what = "character", nlines = 1, sep = "\t", quiet = TRUE)
  colnames(all_data) <- header_row[header_row != ""]

  # d <- data.matrix(t(all_data))
  d <- data.matrix(all_data)
  
  d_filtered <- d[apply(abs(d), 1, max, na.rm = TRUE) > 400000, ]
  d <- d_filtered
  hr <- hclust(dist(d, method="manhattan"), method="complete")
  # hr <- hclust(dist(d, method="euclidean"), method="complete")
  
  out=paste(f,".heatmap.filter400000.png",sep='')
  png(out,res=250,width=2500,height=2500)
  heatmap.2(
    d,
    trace="none",
    #col=c("#FFFFFF", "#000000"),
    col="heat.colors",
    scale="none",
    #RowSideColor=CLASS,
    dendrogram="both",
    cexCol=0.4,
    cexRow=0.4,
    margins=c(15,15)
  )
  dev.off()
  
  out=paste(f,".heatmap.filter400000.svg",sep='')
  svg(out)
  heatmap.2(
    d,
    trace="none",
    #col=c("#FFFFFF", "#000000"),
    scale="none",
    #RowSideColor=CLASS,
    dendrogram="both",
    cexCol=0.4,
    cexRow=0.4,
    margins=c(15,15)
  )
  dev.off()
}




