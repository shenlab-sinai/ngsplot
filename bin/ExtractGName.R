# Extract gene names as ordered by GO algo from ngs.plot output zip file.
# Input should be the file name of the zip file without .zip suffix.

fname <- commandArgs(T)

zip.fname <- paste(fname, 'zip', sep='.')
heatmap.dat <- file.path(fname, 'heatmap.RData')
load(unz(zip.fname, heatmap.dat))

gene.list <- strsplit(go.list[[1]], ':')
gname.list <- sapply(gene.list, function(x) x[1])
write.table(gname.list, file=paste(fname, 'gname.txt', sep='.'), 
            col.names=F, row.names=F)