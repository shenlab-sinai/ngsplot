make_cgi_list <- function(cgi.anno, database, gene.anno=NULL){
	stopifnot(database == 'refseq' || (database == 'ensembl' && !is.null(gene.anno)))
	# Read cpg islands.
	cgi.cols <- c('chrom','start','end','name','score','strand',
		'gname','tid','tstrand','tstart','tend','feature')
	if(database == 'ensembl'){
		cgi.cols[7] <- 'gid'
	}
	cgi <- read.table(cgi.anno, sep="\t", col.names=cgi.cols)
	cgi[cgi$tstrand == "", ]$tstrand <- '+'	# set transcript strand to "+" for intergenic cgi.
	# Format cgi table.
	if(database == 'ensembl'){
		gid_name.tbl <- read.table(gene.anno, sep="\t")
		gid_name.tbl <- unique(gid_name.tbl[, c(4,5)])	# col 4: gid; 5: gname.
		colnames(gid_name.tbl) <- c('gid','gname')
		cgi <- merge(cgi, gid_name.tbl, all.x=T, sort=F)
		cgi.tbl <- data.frame(chrom=cgi$chrom, start=cgi$start, end=cgi$end,
			gid=cgi$gid, gname=cgi$gname, tid=cgi$tid, strand=cgi$tstrand,
			byname.uniq=T, bygid.uniq=T)
	}else{
		cgi.tbl <- data.frame(chrom=cgi$chrom, start=cgi$start, end=cgi$end,
			gid=NA, gname=cgi$gname, tid=cgi$tid, strand=cgi$tstrand,
			byname.uniq=T, bygid.uniq=NA)
	}
	# Split cgi table based on feature annotation.
	split(cgi.tbl, cgi$feature)
}
