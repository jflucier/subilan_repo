```{r}
install.packages("rrcov")
library("rrcov")

if (!("DESeq2" %in% installed.packages())) {
  # Install DESeq2
  BiocManager::install("DESeq2", update = FALSE)
}

if (!requireNamespace('BiocManager', quietly = TRUE))
    install.packages('BiocManager')

BiocManager::install('PCAtools')

install.packages("pheatmap")

```


```{r}

library("rrcov")
library(DESeq2)
library(ggbiplot)
library(devtools)
library("RColorBrewer")
library(corrplot)
library(dplyr)

# # on raw count -- not so good
# f="/storage/Documents/service/externe/ilan/20240308_rnaseq_Khan/DGE/rawCountMatrix.tsv"
# 
# all_data= read.csv(
#   f, 
#   header = TRUE,
#   sep="\t", 
#   row.names=1, 
#   na.strings = "", 
#   stringsAsFactors=FALSE
# )
# 
# head(all_data[-2])
# 
# data <- all_data[, ! names(all_data) %in% "Symbol", drop = F]
# 
# data.cor = cor(data,use = "complete.obs", method = c("spearman"))
# 
# out=paste(f,".all.correlation.png",sep='')
# png(out,res=250,width=2500,height=2500)
# corrplot(
#   data.cor,
#   type = "upper",
#   method = "color",
#   tl.cex = 0.5,
#   title = "ALL group correlation matrix",
#   order = "original",
#   is.corr = T,
#   mar=c(0,0,1,0)
# )
# dev.off()


# on deseq normed counts
f="/storage/Documents/service/externe/ilan/20240308_rnaseq_Khan/DGE/rawCountMatrix.rmoutliers.tsv"
coldata_f="/storage/Documents/service/externe/ilan/20240308_rnaseq_Khan/DGE/coldata.tsv"

all_data= read.csv(
  f, 
  header = TRUE,
  sep="\t", 
  row.names=1, 
  na.strings = "", 
  stringsAsFactors=FALSE
)

head(all_data[-2])

data <- as.matrix(all_data[, ! names(all_data) %in% "Symbol", drop = F])

coldata <- read.csv(
  coldata_f, 
  header = TRUE,
  sep="\t", 
  row.names=1, 
  na.strings = "", 
  stringsAsFactors=FALSE
)
coldata$condition <- factor(coldata$condition)
coldata$repl <- factor(coldata$repl)

head(data)
coldata

all(rownames(coldata) %in% colnames(data))

dds <- DESeqDataSetFromMatrix(data, coldata, design = ~ condition)
featureData <- data.frame(gene=rownames(data))
vsd <- vst(dds, blind=FALSE)
rld <- rlog(dds, blind=FALSE)

head(assay(vsd), 3)


library("RColorBrewer")
library("pheatmap")
dds <- estimateSizeFactors(dds)
# select <- order(rowMeans(counts(dds,normalized=TRUE)),
#                 decreasing=TRUE)[1:20]
# df <- as.data.frame(colData(dds)[,c("condition","repl")])
# pheatmap(assay(vsd)[select,], cluster_rows=FALSE, show_rownames=FALSE,
#          cluster_cols=FALSE, annotation_col=df)

sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$repl, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

out=paste(f,".distmatrix.png",sep='')
png(out,res=250,width=2500,height=2500)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()


library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$repl, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(
  sampleDistMatrix,
  clustering_distance_rows=sampleDists,
  clustering_distance_cols=sampleDists,
  col=colors
)


out=paste(f,".pca.group.png",sep='')
png(out,res=250,width=2500,height=2500)
plotPCA(vsd, intgroup=c("condition"))
dev.off()

out=paste(f,".pca.group.repl.png",sep='')
png(out,res=250,width=2500,height=2500)
plotPCA(vsd, intgroup=c("condition","repl"))
dev.off()



```