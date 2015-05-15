#!/usr/bin/env Rscript

# Extract gene names as ordered by GO algo from ngs.plot output zip file.
# Input should be the file name of the zip file without .zip suffix.

help<-function(){
	cat("Usage: ExtractGname.R file\n")
  cat("Extract gene names and cluster info for ngsplot\n\n")
  cat("File can either be the zip file containing the RData file or the RData 
file directly.\n")
	cat("Output varies depending on input. If input data file has no cluster 
information, only one gene_name.txt file produced for each region.\n")
	cat("If cluster information present, an additional C*.txt file 
produced for each cluster for each region.\n")
	q(status=1)
}

fname <- commandArgs(T)
if(length(fname) != 1) {
	help()
}

require(tools)
fname.ext <- file_ext(fname)
fname.root <- file_path_sans_ext(fname)
if(fname.ext == "zip"){
	heatmap.dat <- file.path(fname.root, 'heatmap.RData')
	load(unz(fname, heatmap.dat))
} else if (fname.ext == "RData") {
  load(fname)
} else { # assume zip file missing .zip suffix. maintains compatibility
  heatmap.dat <- file.path(fname, 'heatmap.RData')
  load(unz(paste(fname, "zip", sep='.'), heatmap.dat))
}

###Step through each gene list
for(i in 1:length(go.list[[1]])){
  split.gname.list <- strsplit(go.list[[1]][[i]], ':')
  gene.tab <- data.frame(do.call(rbind, split.gname.list))
  colnames(gene.tab) <- c("Gene", "Transcript")
  if(!is.na(go.list[[2]][i])) { # with cluster info.
    gene.tab <- data.frame(gene.tab, Cluster=go.list[[2]][[i]])
    gene.cluster.list <- with(gene.tab, split(gene.tab, Cluster, T))
    dumb <- sapply(names(gene.cluster.list), function(cluster) {
      cluster.fname <- paste(fname.root, '.R', i, '.C', cluster, '.csv', 
                             sep='')
      write.csv(gene.cluster.list[[cluster]], row.names=F, file=cluster.fname)
    })
  }
  write.csv(gene.tab, row.names=F,
            file=paste(fname.root, '.R', i, '.gene_name.csv', sep=''))

}
