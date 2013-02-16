#!/usr/bin/env Rscript
#
# Program: ngs.plot.r
# Purpose: Plot sequencing coverages at different genomic regions.
#          Allow overlaying various coverages with gene lists.
# Arguments: coverage file, region to plot, title, output basename.
#            user can also supply a text file describing the interaction between
#            coveage files and gene lists; or a customized region list in BED
#            format.
#
# -- by Li Shen, MSSM
#    Nov 2011.
#

# Deal with command line arguments.
cmd.help <- function(){
	cat("\n")
	cat("Usage: ngs.plot.r -R region_2_plot -C cov_file -O out_base_name [-F reg_further_info] [-D database_choice] [-T title] [-G gene_list] [-I interval_size] [-L flanking_size] [-N flanking_factor]\n")
	cat("       More optional parameters: [-FI forbid_image] [-S random_sample_rate] [-A smooth_function_radius] [-M smooth_method] [-H shaded_area] [-E weigh_genelen] [-P cores_number] [-SE standard_error_points]\n")
	cat("\n## Mandatory parameters:\n")
	cat("  -R     Genomic regions to plot, can be: tss, tes, genebody, exon, cgi or a customized BED file.\n")
	cat("  -C     Coverage file or for multiple plot, a configuration file(must be *.txt, see multiplot.example.txt).\n")
	cat("  -O     Basename for output(WITHOUT suffix). Two files will be generated: pdf image and text file.\n")
	cat("\n## Optional parameters that should be provided in configuration file for multiplot:\n")
	cat("  -G     Gene list to subset regions(default=whole genome).\n")
	cat("  -T     Image title(default=Noname), will be used in figure legend.\n")
	cat("\n## Important optional parameters:\n")
	cat("  -F     If you select genebody, exon or cgi, further information can be provided:\n")
	cat("           for genebody: chipseq(default), rnaseq.\n")
	cat("             RNA-seq is handled differently from ChIP-seq to cope with intronic regions.\n")
	cat("           for exon: canonical(default), variant, promoter, polyA, altAcceptor, altDonor, altBoth.\n")
	cat("             the alt's mean exons with alternative boundaries.\n")
	cat("           for cgi, choose one based on location: ProximalPromoter(default), Promoter1k, Promoter3k,\n")
	cat("             Genebody, Genedesert, OtherIntergenic, Pericentromere.\n")
	cat("  -D     Gene database: refseq(default), ensembl\n")
	cat("  -I     Internal region size(ignored for tss and tes, default for others: genebody=3kb;exon=250;cgi=500;bed=1kb).\n")
	cat("  -L     Flanking region size(default: tss,tes,genebody,bed=1kb;exon=500;cgi=500).\n")
	cat("  -N     Flanking region factor(any numeric >0, will override flanking size if specified).\n")
	cat("           If set, flanking region has floating size with internal size.\n")
	cat("  -S     Randomly sample the regions for plot, must be:(0, 1] (default=1, i.e. all regions).\n")
	cat("  -P     Number of CPUs to be used(default=1). Set 0 to use all detected CPUs.\n")
	cat("  -HM    Number of columns for each heatmap(default=100). Set 0 to turn off heatmap.\n")
	cat("           Generating heatmaps means longer running time and higher memory requirement.\n")
	cat("\n## Misc. parameters:\n")
	cat("  -FI    Forbid image output if set to 1(default=0). This can be useful if you run on a server without graphic output.\n")
	cat("  -M     Smooth method, choose from: mean(default), median\n")
	cat("  -A     Radius used by smooth function, must be:[0, 1)(default=0, i.e. Off). Suggested value: <=0.05.\n")
	cat("           Interpreted as a fraction of the entire plot. Larger value means smoother plot.\n")
	cat("  -H     Use shaded area instead of curves, must be:[0, 1)(default=0, i.e. Off). Suggested value: <0.5.\n")
	cat("           Interpreted as the degree of opacity.\n")
	cat("  -E     Calculate weighted coverage according to internal size, must be: 0 or 1(default=0, i.e. Off).\n")
	cat("           Can be useful, e.g. if you believe longer gene should weigh heavier than shorter gene.\n")
	cat("  -SE    Number of standard errors to calculate. Default=200, set to 0 to turn off.\n")
	cat("\n")
}


###########################################################################
#################### Deal with program input arguments ####################
args <- commandArgs(T)
progpath <- Sys.getenv('NGSPLOT')
if(progpath == ""){
	stop("Set environment variable NGSPLOT before run the program. \
		  See README for details.\n")
}else{
	if(substr(progpath, nchar(progpath), nchar(progpath)) != '/'){	# add trailing slash.
		progpath <- paste(progpath, '/', sep='')
	}
}
source(paste(progpath, 'lib/parse.args.r', sep=''))
args.tbl <- parse.args(args, c('-C', '-R', '-O'))
if(is.null(args.tbl)){
	cmd.help()
	stop('Error in parsing command line arguments. Stop.\n')
}
covfile <- args.tbl['-C']
reg2plot <- args.tbl['-R']
basename <- args.tbl['-O']

# Determine coverage-title-genelist relationship.
if(length(grep('.txt$', covfile))>0){
	ctg.tbl <- read.table(covfile, sep="\t", 
				col.names=c('cov','glist','title'), as.is=T)
}else{
	if('-G' %in% names(args.tbl)){
		glist <- args.tbl['-G']
	}else{
		glist <- '-1'
	}
	if('-T' %in% names(args.tbl)){
		title <- args.tbl['-T']
	}else{
		title <- 'Noname'
	}
	ctg.tbl <- data.frame(cov=covfile, glist=glist, title=title, 
							stringsAsFactors=F)
}

# Collapse coverage files to speed up loading.
cov.u <- unique(ctg.tbl$cov)


# Load the 1st coverage. The genome name is then used for loading gene models.
cov2load <- cov.u[1]
load(cov2load)	# genome, nreads, read.coveage, read.coveage.n
if(genome == 'mm9'){	# load genome: refseq, ensembl, cgiUCSC
	load(paste(progpath, 'database/mm9.RData', sep=''))
}else if(genome == 'hg19'){
	load(paste(progpath, 'database/hg19.RData', sep=''))
}else if(genome == 'rn4'){
	load(paste(progpath, 'database/rn4.RData', sep=''))
}else{
	stop(paste('Unsupported genome: ', genome, '. Stop.\n', sep=''))
}

#### Image output forbidden tag. ####
if('-FI' %in% names(args.tbl)){	
	stopifnot(as.integer(args.tbl['-FI']) >= 0)
	fi_tag <- as.integer(args.tbl['-FI'])
}else{
	fi_tag <- as.integer(0)
}

#### Create heatmap tag. ####
if('-HM' %in% names(args.tbl)){	
	stopifnot(as.integer(args.tbl['-HM']) >= 0)
	hm_cols <- as.integer(args.tbl['-HM'])
}else{
	hm_cols <- as.integer(100)
}

#### Database. ####
if('-D' %in% names(args.tbl)){	
	database <- as.character(args.tbl['-D'])
	if(database == 'refseq'){	# choose database.
		genemodel <- refseq
	}else if(database == 'ensembl'){
		genemodel <- ensembl
	}else{
		stop('Unsupported database.')
	}
}else{
	genemodel <- refseq
}

# Function to remove ".fa" if needed, or to add "chr" if absent.
chromFormat <- function(crn, ...){
	crn <- sub('.fa$', '', crn)
	nochr.i <- grep('^chr', crn, invert=T)
	crn[nochr.i] <- paste('chr', crn[nochr.i], sep='')
	crn
}

#### Determine the set of genomic coordinates. ####
if(reg2plot == 'tss' || reg2plot == 'tes'){	
	genome.coord <- genemodel$genebody
}else if(reg2plot == 'genebody'){
	if('-F' %in% names(args.tbl)){
		finfo <- as.character(args.tbl['-F'])
		gb.allowed <- c('chipseq', 'rnaseq')
		stopifnot(finfo %in% gb.allowed)
	}else{
		finfo <- 'chipseq'
	}
	genome.coord <- genemodel$genebody
}else if(reg2plot == 'exon'){
	if('-F' %in% names(args.tbl)){
		finfo <- as.character(args.tbl['-F'])
		exon.allowed <- c('canonical', 'variant', 'promoter', 'polyA', 
							'altAcceptor', 'altDonor', 'altBoth')
		stopifnot(finfo %in% exon.allowed)
	}else{
		finfo <- 'canonical'
	}
	genome.coord <- genemodel$exon
}else if(reg2plot == 'cgi'){
	if('-F' %in% names(args.tbl)){
		finfo <- as.character(args.tbl['-F'])
		cgi.allowed <- c("Genebody", "Genedesert", "OtherIntergenic", 
							"Pericentromere", "Promoter1k", "Promoter3k", 
							"ProximalPromoter")
		stopifnot(finfo %in% cgi.allowed)
	}else{
		finfo <- 'ProximalPromoter'
	}
	genome.coord <- genemodel$cgi
}else{	# if not above, then BED file.
	bed.coord <- read.table(reg2plot, sep="\t")
	if(ncol(bed.coord) <3){
		stop('Input file must contain at least 3 columns!')
	}
	genome.coord <- data.frame(chrom=chromFormat(bed.coord[, 1]), 
						start=bed.coord[, 2]+1, end=bed.coord[, 3], 
						gid=NA, gname='N', tid='N', strand='+', 
						byname.uniq=T, bygid.uniq=NA)
	if(ncol(bed.coord) >=4){
		genome.coord$gname <- bed.coord[, 4]
	}
	if(ncol(bed.coord) >=5){
		genome.coord$tid <- bed.coord[, 5]
	}
	if(ncol(bed.coord) >=6){
		genome.coord$strand <- bed.coord[, 6]
	}
	reg2plot <- 'bed'	# rename for information retrieval
}
if(reg2plot == 'genebody' && finfo == 'rnaseq'){
	rnaseq.gb <- T
}else{
	rnaseq.gb <- F
}
if(reg2plot == 'exon' || reg2plot == 'cgi'){	# subset specific region.
	genome.coord <- genome.coord[[finfo]]
}

#### Interval region size. ####
if('-I' %in% names(args.tbl)){	
	stopifnot(as.integer(args.tbl['-I']) > 0)
	intsize <- args.tbl['-I']
}else{
	int.tbl <- c(3000,250,500,1000)
	names(int.tbl) <- c('genebody','exon','cgi','bed')
	intsize <- int.tbl[reg2plot]
        if(reg2plot == 'bed' && genome.coord$end == genome.coord$start){
	        intsize <- 1
        }
}
if(reg2plot == 'tss' || reg2plot == 'tes'){
	intsize <- 1
}
intsize <- as.integer(intsize)

#### Flanking region size. ####
if('-L' %in% names(args.tbl)){	
	stopifnot(as.integer(args.tbl['-L']) >= 0)
	flanksize <- args.tbl['-L']
}else{
	flank.tbl <- c(1000,1000,1000,500,500,1000)
	names(flank.tbl) <- c('tss','tes','genebody','exon','cgi','bed')
	flanksize <- flank.tbl[reg2plot]
}
flanksize <- as.integer(flanksize)

#### Flanking size factor. ####
if('-N' %in% names(args.tbl) && !('-L' %in% names(args.tbl))){	
	stopifnot(as.numeric(args.tbl['-N']) >= 0)
	flankfactor <- as.numeric(args.tbl['-N'])
	flanksize <- floor(intsize*flankfactor)
}else{
	flankfactor <- 0.0
}
if(rnaseq.gb){	# RNA-seq plotting.
	flanksize <- 0
}

#### Random sampling rate. ####
if('-S' %in% names(args.tbl)){	
	samprate <- as.numeric(args.tbl['-S'])
	stopifnot(samprate > 0 && samprate <= 1)
	recs <- which(genome.coord$byname.uniq)	# records to sample from.
	if(samprate < 1){	# random sample indices.
		samp.i <- sample(recs, floor(samprate*length(recs)))
	}
}else{
	samprate <- 1.0
}

#### Shaded area alpha. ####
if('-H' %in% names(args.tbl)){	
	shade.alp <- as.numeric(args.tbl['-H'])
	stopifnot(shade.alp >= 0 && shade.alp < 1)
}else{
	shade.alp <- 0
}

#### Smooth function radius. ####
if('-A' %in% names(args.tbl)){	
	smooth.radius <- as.numeric(args.tbl['-A'])
	stopifnot(smooth.radius >= 0 && smooth.radius < 1)
}else{
	smooth.radius <- .0
}

#### Smoothing method. ####
if('-M' %in% names(args.tbl)){	
	smooth.method <- as.character(args.tbl['-M'])
	stopifnot(smooth.method == 'mean' || smooth.method == 'median')
}else{
	smooth.method <- 'mean'
}

#### Weighted coverage. ####
if('-E' %in% names(args.tbl)){	
	stopifnot(as.integer(args.tbl['-E']) >= 0)
	weight.genlen <- as.integer(args.tbl['-E'])
}else{
	weight.genlen <- as.integer(0)
}

##### Set cores number. ####
if('-P' %in% names(args.tbl)){
	stopifnot(as.integer(args.tbl['-P']) >= 0)
	cores.number <- as.integer(args.tbl['-P'])
}else{
	cores.number <- as.integer(1)
}

#### Set number of points to calculate standard errors. ####
if('-SE' %in% names(args.tbl)){	
	stopifnot(as.integer(args.tbl['-SE']) >= 0)
	stderror.number <- as.integer(args.tbl['-SE'])
	confiMat <- matrix(0, nrow=stderror.number+1, ncol=nrow(ctg.tbl))
}else{
	stderror.number <- as.integer(200)
	confiMat <- matrix(0, nrow=stderror.number+1, ncol=nrow(ctg.tbl))
}
colnames(confiMat) <- ctg.tbl$title
if(stderror.number == 0){
	confiMat <- NULL
}

# Create the matrix to store plotting data.
if(rnaseq.gb){
	regcovMat <- matrix(0, nrow=intsize, ncol=nrow(ctg.tbl))
}else{
	regcovMat <- matrix(0, nrow=2*flanksize + intsize, ncol=nrow(ctg.tbl))
}
colnames(regcovMat) <- ctg.tbl$title

############### End arguments configuration #####################
#################################################################




# Load required libraries.
require(ShortRead) || {
	source("http://bioconductor.org/biocLite.R")
	biocLite(ShortRead)
	require(ShortRead)
}
require(BSgenome) || {
	source("http://bioconductor.org/biocLite.R")
	biocLite(BSgenome)
	require(BSgenome)
}
require(doMC)

if(cores.number == 0){
	registerDoMC()
} else {
	registerDoMC(cores.number)
}


source(paste(progpath, 'lib/coverage.r', sep=''))
source(paste(progpath, 'lib/plotlib.r', sep=''))


# Setup SEM sample points.
if(stderror.number > 0){
	stderror.pos <- round(seq(1, nrow(regcovMat), length.out=stderror.number+1))
}

# Function to calculate standard error
calcSem <- function(x){ sd(x, na.rm=T)/sqrt(length(x)) }


# Genomic enrichment for all profiles in the config. Use this for heatmaps.
# Pre-allocate list size to improve efficiency.
enrichList <- as.list(rep(NA, nrow(ctg.tbl)))


old_flanksize <- flanksize # still needed?

# Go through all remaining unique coverage files.
# for(i in 2:length(cov.u)){
i <- 1
while(i <= length(cov.u)){
	# "same.cov.r" contains the row numbers of the config file.
	same.cov.r <- which(ctg.tbl$cov == cov2load) # "cov2load" is set before. 

	# Go through all gene lists associated with each coverage.
	# "r" is also the column position of the plot matrix.
  	foreach(r=iter(same.cov.r)) %do% {

		lname <- ctg.tbl$glist[r]	# retrieve gene list names.
		if(lname == '-1'){	# use genome as gene list.
			if(samprate < 1){
				plot.coord <- genome.coord[samp.i, ]
			}else{
				plot.coord <- subset(genome.coord, byname.uniq)
			}
		}else{	# read gene list from text file.
			gene.list <- read.table(lname, as.is=T, comment.char='#')$V1
			subset.idx <- c(which(genome.coord$gname %in% gene.list & 
									genome.coord$byname.uniq),
							which(genome.coord$tid %in% gene.list)
							)
			if(!all(is.na(genome.coord$gid))){
				subset.idx <- c(subset.idx, 
								which(genome.coord$gid %in% gene.list & 
									genome.coord$bygid.uniq)
								)
			}
			plot.coord <- genome.coord[subset.idx, ]
		}

		# Extract coverage and combine into a matrix.
		# result.matrix <- foreach(k=1:nrow(plot.coord), 
		# 				.combine='rbind', .multicombine=F) %dopar% {
		cov.result <- foreach(k=1:nrow(plot.coord)) %dopar% {
        	do.par.cov(k, plot.coord, read.coverage.n, rnaseq.gb,
				flankfactor, reg2plot, genemodel, weight.genlen,
				intsize, old_flanksize, flanksize)
		}
		result.matrix <- do.call('rbind', cov.result)

		# Calc avg. profile.
		regcovMat[, r] <- apply(result.matrix, 2, function(x) mean(x, na.rm=T))	

		# Calculate SEM if needed. Shut off SEM in single gene case.
		if(nrow(result.matrix) > 1 && stderror.number > 0){
			confiMat[, r] <- apply(result.matrix[, stderror.pos], 2, calcSem)
		}

		# Sample and book-keep this matrix for heatmap.
		if(hm_cols){
			enrichList[[r]] <- spline_mat(result.matrix, hm_cols)
		}
	}
	# load a new coverage file.
	i <- i + 1
	if(i <= length(cov.u)){
		cov2load <- cov.u[i]
		load(cov2load)
	}
}
flanksize <- old_flanksize	# recover the original flanksize.

# Smooth plot if specified.
if(smooth.radius > 0){
	regcovMat <- smoothplot(regcovMat, smooth.radius, smooth.method)
}

# Create image file and plot data into it.
if(!fi_tag){
	# Average profile plot.
	default.width <- 8	# in inches.
	default.height <- 7
	out.plot <- paste(basename, '.pdf', sep='')
	pdf(out.plot, width=default.width, height=default.height)
	xticks <- gen_xticks(reg2plot, intsize, flanksize, flankfactor)
	plotmat(regcovMat, ctg.tbl$title, xticks, shade.alp, confiMat)
	dev.off()

	# Heatmap.
	if(hm_cols){
		# Identify unique regions.
		reg.list <- as.factor(ctg.tbl$glist)
		uniq.reg <- levels(reg.list)
		# Number of plots per region.
		reg.np <- sapply(uniq.reg, function(r) sum(reg.list==r))
		# Number of genes per region.
		reg.ng <- sapply(uniq.reg, function(r){
			ri <- which(reg.list==r)[1]
			nrow(enrichList[[ri]])
		})
		# Setup image size.
		unit.width <- 4	# in inches.
		reduce.ratio <- 20 	# col to row reduction because #gene is large.
		hm.width <- unit.width * max(reg.np)
		ipl <- .2 # inches per line. Obtained from par->'mai', 'mar'.
		m.bot <- 2; m.lef <- 1; m.top <- 2; m.rig <- 1 # margin size in lines.
		# Convert #gene to image height.
		reg.hei <- sapply(reg.ng, function(r){
			r * unit.width / hm_cols / reduce.ratio + 
				m.bot * ipl + m.top * ipl # margins are like intercepts.
		})
		hm.height <- sum(reg.hei)

		# Setup output device.
		out.hm <- paste(basename, '.hm.pdf', sep='')
		pdf(out.hm, width=hm.width, height=hm.height)
		par(mar=c(m.bot, m.lef, m.top, m.rig))

		# Setup layout of the heatmaps.
		lay.mat <- matrix(0, ncol=max(reg.np), nrow=length(reg.np))
		f.sta <- 1 # figure number start.
		for(i in 1:length(reg.np)){
			f.end <- f.sta + reg.np[i] - 1 # figure number end.
			lay.mat[i, 1:reg.np[i]] <- f.sta : f.end
			f.sta <- f.sta + reg.np[i]
		}
		layout(lay.mat, heights=reg.hei)

		plotheat(reg.list, uniq.reg, enrichList, ctg.tbl$title, xticks)
		dev.off()
	}
}

# Save plotting data of reads density to a text file.
out.txt <- paste(basename, '.txt', sep='')
out.header <- c(
	'#Do NOT change the following lines if you want to re-draw the image with \
	 replot.r! If you change the matrix values, cut and paste the commented \
	 lines into your new file and run replot.r.',
	paste('#reg2plot:', reg2plot, sep=''),
	paste('#flanksize:', flanksize, sep=''),
	paste('#intsize:', intsize, sep=''),
	paste('#flankfactor:', flankfactor, sep=''),
	paste('#shade.alp:', shade.alp, sep=''),
	paste('#rnaseq.gb:', rnaseq.gb, sep=''),
	paste('#width:', default.width, sep=''),
	paste('#height:', default.height, sep=''))
writeLines(out.header, out.txt)
suppressWarnings(write.table(regcovMat, append=T, file=out.txt, 
					row.names=F, sep="\t", quote=F))

# Save plotting data of standard error to a text file.
if(!is.null(confiMat)){
	out.txt <- paste(basename, '_stderror.txt', sep='')
	out.header <- '#Do NOT change the following lines if you want to re-draw \
					the image with replot.r!'
	writeLines(out.header, out.txt)
	suppressWarnings(write.table(confiMat, append=T, file=out.txt, 
						row.names=F, sep="\t", quote=F))
}
