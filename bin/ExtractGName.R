#!/usr/bin/env Rscript

# Extract gene names as ordered by GO algo from ngs.plot output zip file.
# Input should be the file name of the zip file without .zip suffix.
library(stringr)

fname <- commandArgs(T)

zip.fname <- paste(fname, 'zip', sep='.')
heatmap.dat <- file.path(fname, 'heatmap.RData')
load(unz(zip.fname, heatmap.dat))

gene.list <- as.data.frame(go.list, stringsAsFactors = FALSE)
gene.list<-as.data.frame(list(str_split_fixed(gene.list[,1],":",2),gene.list[,2]))

write.table(gene.list, file=paste(fname, 'gname.txt', sep='.'), 
            col.names=F, row.names=F)

clusters <- max(gene.list[,3])
for( i in 1:clusters ){
  filename <- paste("cluster",i,".txt",sep="")
  write.table(gene.list[gene.list[,3]==i,1], file=paste(fname, filename, sep='.'),  col.names=F, row.names=F, quote = FALSE)
}

