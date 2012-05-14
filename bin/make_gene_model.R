make_gene_model <- function(txtfile){
	stopifnot(is.character(txtfile))
	gene.tbl <- read.table(txtfile, comment.char='', sep="\t", as.is=T, col.names=c('chrom', 'start', 'end', 'gid', 'gname', 'tid', 'strand', 'region', 'class'))
	# genebody.
	genebody <- subset(gene.tbl, region=='genebody')
	# choose longest transcript for each gene name.
	gb.byname.sorted <- genebody[order(genebody$gname, (genebody$end - genebody$start), decreasing=T), ]
	byname.uniq.tid <- gb.byname.sorted[!duplicated(gb.byname.sorted$gname), ]$tid
	genebody <- cbind(genebody[, 1:7], byname.uniq=(genebody$tid %in% byname.uniq.tid))
	# choose longest transcript for each gene id. 
	if(all(is.na(genebody$gid))){
		bygid.uniq.tid <- NULL
		genebody <- cbind(genebody, bygid.uniq=NA)
	}else{
		gb.bygid.sorted <- genebody[order(genebody$gid, (genebody$end - genebody$start), decreasing=T), ]
		bygid.uniq.tid <- gb.bygid.sorted[!duplicated(gb.bygid.sorted$gid), ]$tid
		genebody <- cbind(genebody, bygid.uniq=(genebody$tid %in% bygid.uniq.tid))
	}
	
	# exon.
	exon <- subset(gene.tbl, region=='exon')
	exon <- cbind(exon[, -8], byname.uniq=(exon$tid %in% byname.uniq.tid))
	if(all(is.na(exon$gid))){
		exon <- cbind(exon, bygid.uniq=NA)
	}else{
		exon <- cbind(exon, bygid.uniq=(exon$tid %in% bygid.uniq.tid))
	}
	exon.class <- split(exon, exon$class)

	# exonmodel.
	exon.tid <- split(exon, exon$tid)
	library(IRanges)
	exonmodel <- lapply(exon.tid, function(tid) {
		tid <- tid[order(tid$start), ]
		list(ranges=IRanges(start=tid$start, end=tid$end))
	})
	exonmodel <- exonmodel[match(genebody$tid, names(exonmodel))]
	
	# return a nested list.
	list(genebody=genebody,
		exon=exon.class,
		exonmodel=exonmodel)
}

# Example lines from Ensembl GTF.
#chr6	108733053	108773717	ENSMUSG00000030105	Arl8b	ENSMUST00000032196	+	genebody	NA
#chr6	108733053	108733371	ENSMUSG00000030105	Arl8b	ENSMUST00000032196	+	exon	promoter
#chr6	108763621	108763701	ENSMUSG00000030105	Arl8b	ENSMUST00000032196	+	exon	canonical
#chr6	108764978	108765051	ENSMUSG00000030105	Arl8b	ENSMUST00000032196	+	exon	canonical
#chr6	108765252	108765345	ENSMUSG00000030105	Arl8b	ENSMUST00000032196	+	exon	canonical
#chr6	108768239	108768306	ENSMUSG00000030105	Arl8b	ENSMUST00000032196	+	exon	canonical
#chr6	108768537	108768607	ENSMUSG00000030105	Arl8b	ENSMUST00000032196	+	exon	canonical
#chr6	108771497	108773717	ENSMUSG00000030105	Arl8b	ENSMUST00000032196	+	exon	polyA
#chr3	92054654	92093169	ENSMUSG00000074445	Sprr2a3	ENSMUST00000090872	+	genebody	NA
#chr3	92054654	92054696	ENSMUSG00000074445	Sprr2a3	ENSMUST00000090872	+	exon	promoter
