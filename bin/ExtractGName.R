#!/usr/bin/env Rscript

# Extract gene names as ordered by GO algo from ngs.plot output zip file.
# Input should be the file name of the zip file without .zip suffix.
library(stringr)

fname <- commandArgs(T)

zip.fname <- paste(fname, 'zip', sep='.')
heatmap.dat <- file.path(fname, 'heatmap.RData')
load(unz(zip.fname, heatmap.dat))

gene.list <- as.data.frame(go.list, stringsAsFactors = FALSE)
gname.list <- str_split_fixed(gene.list[,1],':',2)[,1]
write.table(gname.list, file=paste(fname, 'gname.txt', sep='.'), 
            col.names=F, row.names=F)
