#!/usr/bin/env Rscript

# Extract gene names as ordered by GO algo from ngs.plot output zip file.
# Input should be the file name of the zip file without .zip suffix.

help<-function(){
	cat("Usage: ExtractGname.R file\nExtract function for ngsplot\nFile can either be the zip file containing the RData file or the RData file directly\n")
	q(status=1)
}

fname <- commandArgs(T)
fname<-"tmp.RData"
if( length(fname) < 1 || length(fname) > 1){
	help()
}

if(grepl(".zip",fname)){
  root<-strsplit(fname,'.',T)[[1]][1]
	heatmap.dat <- file.path(root, 'heatmap.RData')
	load(unz(fname, heatmap.dat))
} else if (grepl(".RData",fname)) {
  root<-strsplit(fname,'.',T)[[1]][1]
  load(fname)
} else { #assume zip file missing .zip suffix. maintains compatibility
  heatmap.dat <- file.path(fname, 'heatmap.RData')
  load(unz(paste(fname,"zip",sep='.'), heatmap.dat))
}

###Step through each gene list
for(i in 1:length(go.list[[1]])){
  if(is.na(go.list[[2]][i])){ #no cluster information
    gene.list <- data.frame(go.list[[1]][i], stringsAsFactors = FALSE)
    
    gname.list <- strsplit(gene.list[,1], ':')
    gene.list <- data.frame(do.call(rbind, gname.list))

    write.table(gene.list, file=paste(root, i,'gene_name.txt', sep='.'),
                col.names=F, row.names=F)
  } else{ # With Clusters
    gene.list <- data.frame(gene=go.list[[1]][i],cluster=go.list[[2]][i], stringsAsFactors = FALSE)
    
    gname.list <- strsplit(gene.list[,1], ':')
    gene.list <- data.frame(do.call(rbind, gname.list),gene.list[,2])
    
    write.table(gene.list, file=paste(root, i,'gene_name.txt', sep='.'),
                col.names=F, row.names=F)
    clusters <- max(gene.list[,3])
    for( j in 1:clusters ){ # Write Clusters
      filename <- paste("cluster",j,".txt",sep="")
      write.table(gene.list[gene.list[,3]==j,1], file=paste(root,i, filename, sep='.'),  col.names=F, row.names=F, quote = FALSE)
    }
  }
}
