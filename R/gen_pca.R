#!/usr/bin/env Rscript

library("optparse")

option_list = list(
  make_option(c("-gm", "--gene_matrix"), type="character", default=NULL,
              help="diann gene matrix output file name", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$gene_matrix)){
  print_help(opt_parser)
  stop("script requires a diann gene matrix file", call.=FALSE)
}

f <- opt$gene_matrix
# f <- "/storage/Documents/service/externe/ilan/20240702_mouse_ms_organoid/report.unique_genes_matrix.reformat.tsv"

message ("Preparing data")
all_data= read.csv(
  f,
  header = TRUE,
  sep="\t",
  row.names=1,
  na.strings = "",
  stringsAsFactors=FALSE
)
# num_data <- data.frame(data.matrix(all_data))
# numeric_columns <- sapply(num_data,function(x){mean(as.numeric(is.na(x)))<0.5})
# final_data <- data.matrix(num_data[,numeric_columns], all_data[,!numeric_columns])
filtered_data = all_data[,c(2:ncol(all_data))]
my_group <- all_data[,1]

message ("Generate PCA")

library(pcaMethods)

# PCA with NA solution taken from https://stackoverflow.com/questions/49641896/r-how-to-use-ggbiplot-with-pcares-object-plot-pca-results-of-data-with-missing
filtered_data_noNA = filtered_data
filtered_data_noNA[is.na(filtered_data_noNA)]<-0
pca.obj <- prcomp(filtered_data_noNA, center=TRUE, scale.=TRUE)
pca.obj2 <- pca(filtered_data, method="svdImpute", center=TRUE, scale.=TRUE)

# pca.obj$x<-pca.obj2@scores 
# pca.obj$rotation<-pca.obj2@loadings 
# pca.obj$sdev<-pca.obj2@sDev
# pca.obj$center<-pca.obj2@center
# pca.obj$scale<-pca.obj2@scale

library("ggplotify")
library("ggbiplot")

filtered_my_group = all_data[,1]
cols <- c("KO" = "blue", "WT" = "red")
x = as.factor(filtered_my_group)

out=paste(f,".pca1-2.png",sep='')
png(out,res=250,width=2500,height=2500)
# out=paste(f,".pca1-2.svg",sep='')
# svg(out)
print(
  ggbiplot(
    pca.obj,
    choices = 1:2,
    groups=x,
    var.axes=F,
    ellipse=T,
    labels.size = 24,
    varname.size = 24
  )
)+
  scale_colour_manual(values = cols) +
  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

out=paste(f,".pca1-3.png",sep='')
png(out,res=250,width=2500,height=2500)
print(
  ggbiplot(
    pca.obj,
    choices = c(1,3),
    groups=x,
    var.axes=F,
    ellipse=T,
    labels.size = 24,
    varname.size = 24
  )
)+
  scale_colour_manual(values = cols) +
  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()

out=paste(f,".pca2-3.png",sep='')
png(out,res=250,width=2500,height=2500)
print(
  ggbiplot(
    pca.obj,
    choices = c(2,3),
    groups=x,
    var.axes=F,
    ellipse=F,
    labels.size = 24,
    varname.size = 24
  )
)+
  scale_colour_manual(values = cols) +
  theme(legend.direction = 'horizontal', legend.position = 'top')
dev.off()




