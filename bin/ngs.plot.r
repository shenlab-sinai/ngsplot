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
# Created:	  	Nov 2011.
# Last Updated: Feb 2013.
#


# Deal with command line arguments.
cmd.help <- function(){
	cat("\n")
	cat("Usage: ngs.plot.r -R region -C [cov|config]_file -O out_basename \
		[Optional]\n")
	cat("\n## Mandatory parameters:\n")
	cat("  -R   Genomic regions to plot: tss, tes, genebody, exon, cgi or *.bed\n")
	cat("  -C   Coverage file or a configuration file for multiplot\n")
	cat("  -O   Basename for output\n")
	cat("## Optional parameters related to configuration file:\n")
	cat("  -G   Gene list to subset regions\n")
	cat("  -T   Image title\n")
	cat("## Important optional parameters:\n")
	cat("  -F   Further information can be provided to subset regions:\n")
	cat("         for genebody: chipseq(default), rnaseq.\n")
	cat("         for exon: canonical(default), variant, promoter, polyA,\n")
	cat("           altAcceptor, altDonor, altBoth.\n")
	cat("         for cgi: ProximalPromoter(default), Promoter1k, Promoter3k,\n")
	cat("           Genebody, Genedesert, OtherIntergenic, Pericentromere.\n")
	cat("  -D   Gene database: refseq(default), ensembl\n")
	cat("  -I   Shall interval be larger than flanking in plot?(0 or 1)\n")
	cat("  -L   Flanking region size\n")
	cat("  -N   Flanking region factor(will override flanking size)\n")
	cat("  -S   Randomly sample the regions for plot, must be:(0, 1]\n")
	cat("  -P   #CPUs to be used. Set 0 to use all detected\n")
	# cat("  -PT  #data points to be used in plot(default=100)\n")
	cat("## Misc. parameters:\n")
	cat("  -GO  Gene order algorithm used in heatmaps: hc(default), total, max,\n")
	cat("         prod, diff, pca and none(according to gene list supplied)\n")
	cat("  -SE  Shall standard errors be plotted?(0 or 1)\n")
	cat("  -RB  The fraction of extreme values to be trimmed on both ends\n")
	cat("         default=0.05, set 0 to keep all data\n")
	cat("  -FC  Flooding fraction:[0, 1), default=0.02\n")
	cat("  -FI  Forbid image output if set to 1(default=0)\n")
	cat("  -A   Radius used by smooth function, suggested value:[0, 0.05]\n")
	cat("         default=0, i.e. no smoothing.")
	cat("  -M   Smooth method: mean(default) or median\n")
	cat("  -H   Opacity of shaded area, suggested value:[0, 0.5]\n")
	cat("         default=0, i.e. no shading, just curves\n")
	# cat("  -E   Calculate weighted coverage according to region size(experimental!)\n")
	cat("\n")
}


###########################################################################
#################### Deal with program input arguments ####################
args <- commandArgs(T)

# Program environment variable.
progpath <- Sys.getenv('NGSPLOT')
if(progpath == ""){
	stop("Set environment variable NGSPLOT before run the program. \
		  See README for details.\n")
}else{  
	if(substr(progpath, nchar(progpath), nchar(progpath)) != '/'){	
		progpath <- paste(progpath, '/', sep='')  # add trailing slash.
	}
}

# Input argument parser.
source(paste(progpath, 'lib/parse.args.r', sep=''))
args.tbl <- parse.args(args, c('-C', '-R', '-O'))
if(is.null(args.tbl)){
	cmd.help()
	stop('Error in parsing command line arguments. Stop.\n')
}

# Load required libraries.
suppressMessages(require(ShortRead, warn.conflicts=F)) || {
	source("http://bioconductor.org/biocLite.R")
	biocLite(ShortRead)
	if(!suppressMessages(require(ShortRead, warn.conflicts=F))){
		stop('Loading package ShortRead failed!')
	}
}
suppressMessages(require(BSgenome, warn.conflicts=F)) || {
	source("http://bioconductor.org/biocLite.R")
	biocLite(BSgenome)
	if(!suppressMessages(require(BSgenome, warn.conflicts=F))){
		stop('Loading package BSgenome failed!')
	}
}
suppressMessages(require(doMC, warn.conflicts=F)) || {
	install.packages('doMC')
	if(!suppressMessages(require(doMC, warn.conflicts=F))){
		stop('Loading package doMC failed!')
	}
}

# Configuration: coverage-genelist-title relationships.
ctg.tbl <- ConfigTbl(args.tbl)

# Collapse coverage files to speed up loading.
cov.u <- unique(ctg.tbl$cov)

# Load the 1st coverage. The genome name is then used for loading gene models.
cov2load <- cov.u[1]
load(cov2load)
# if(genome == 'mm9'){
# 	load(paste(progpath, 'database/mm9.RData', sep=''))
# }else if(genome == 'hg19'){
# 	load(paste(progpath, 'database/hg19.RData', sep=''))
# }else if(genome == 'rn4'){
# 	load(paste(progpath, 'database/rn4.RData', sep=''))
# }else{
# 	stop(paste('Unsupported genome: ', genome, '. Stop.\n', sep=''))
# }


# Setup variables from arguments.
argvar.list <- setup_vars(args.tbl, ctg.tbl)
attach(argvar.list)
# This sets a large number of variables for use by the following procedures.
# reg2plot: tss, tes, genebody, *.bed...
# basename: output file basename
# fi_tag: tag for forbidding image output
# lgint: tag for large interval
# flanksize: flanking region size
# flankfactor: flanking region factor
# samprate: sampling rate
# shade.alp: shade area alpha
# smooth.radius: smooth function radius
# smooth.method: as is
# cores.number: #CPUs
# se: tag for plotting stand errors
# robust: robust stat fraction
# flood.frac: flooding fraction.
# go.algo: gene order algorithm used in heatmaps.


# Register doMC with CPU number.
if(cores.number == 0){
	registerDoMC()
} else {
	registerDoMC(cores.number)
}


# Setup plot-related coordinates and variables.
source(paste(progpath, 'lib/genedb.r', sep=''))
plotvar.list <- SetupPlotCoord(args.tbl, ctg.tbl, progpath, genome, reg2plot, 
								samprate)
attach(plotvar.list)
# This sets plot coordinates and a few variables.
# coord.list: a list of coordinates as data.frame for all unique regions.
# rnaseq.gb: tag for RNA-seq data.
# reg.list: region list as factor vector.
# uniq.reg: unique region list as character list.
# pint: tag for point interval.
# exonmodel: exon ranges if rnaseq.gb=True.


# Setup data points for plot.
source(paste(progpath, 'lib/plotlib.r', sep=''))
pts.list <- SetPtsSpline(pint)
attach(pts.list)
# pts: data points for avg. profile and standard errors.
# m.pts: middle data points. For pint, m.pts=1.
# f.pts: flanking region data points.

# Setup matrix for avg. profiles.
regcovMat <- CreatePlotMat(pts, ctg.tbl)

# Setup matrix for standard errors.
confiMat <- CreateConfiMat(se, pts, ctg.tbl)

# Genomic enrichment for all profiles in the config. Use this for heatmaps.
# Pre-allocate list size to improve efficiency.
enrichList <- as.list(rep(NA, nrow(ctg.tbl)))

###################################################
# Here start to extract coverages for all genomic regions and calculate 
# data for plotting.

# Load coverage extraction lib.
source(paste(progpath, 'lib/coverage.r', sep=''))

i <- 1  # cov2load was previously loaded. 
# Go through all remaining unique coverage files.
while(i <= length(cov.u)){
	# "same.cov.r" contains the row numbers of the config file.
	same.cov.r <- which(ctg.tbl$cov == cov2load) # "cov2load" is set before. 

	# Go through all gene lists associated with each coverage.
	# "r" is also the column position of the plot matrix.
  	foreach(r=iter(same.cov.r)) %do% {
		reg <- ctg.tbl$glist[r]	# retrieve gene list names.

		# Extract coverage and combine into a matrix.
		result.matrix <- foreach(rec=iter(coord.list[[reg]], by='row'), 
							.combine='rbind', .multicombine=T, 
							.maxcombine=1000) %dopar% {
			# Convert coord record to exon ranges.
			if(rnaseq.gb) {
				exon.ranges <- exonmodel[[rec$tid]]$ranges
			} else {
				exon.ranges <- NULL
			}
        	do.par.cov(rec, read.coverage.n, flanksize, flankfactor, 
        				m.pts, f.pts, reg2plot, pint, exon.ranges)
		}
		# result.matrix <- do.call('rbind', cov.result)

		# Calc avg. profile.
		regcovMat[, r] <- apply(result.matrix, 2, function(x) mean(x, 
								trim=robust, na.rm=T))	

		# Calculate SEM if needed. Shut off SEM in single gene case.
		if(nrow(result.matrix) > 1 && se){
			confiMat[, r] <- apply(result.matrix, 2, function(x) {
									CalcSem(x, robust)
								})
		}

		# Sample and book-keep this matrix for heatmap.
		enrichList[[r]] <- result.matrix
	}
	# load a new coverage file.
	i <- i + 1
	if(i <= length(cov.u)){
		cov2load <- cov.u[i]
		load(cov2load)
	}
}

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
	xticks <- gen_xticks(reg2plot, pint, lgint, pts, flanksize, flankfactor)

	plotmat(regcovMat, ctg.tbl$title, xticks, pts, m.pts, f.pts, pint,
			shade.alp, confiMat)
	dev.off()

	# Heatmap.
	# Identify unique regions.
	# reg.list <- as.factor(ctg.tbl$glist)
	# uniq.reg <- levels(reg.list)
	# Number of plots per region.
	reg.np <- sapply(uniq.reg, function(r) sum(reg.list==r))
	# Number of genes per region.
	reg.ng <- sapply(uniq.reg, function(r){
					ri <- which(reg.list==r)[1]
					nrow(enrichList[[ri]])
				})
	# Setup image size.
	unit.width <- 4	# in inches.
	reduce.ratio <- 10 	# col to row reduction because #gene is large.
	hm.width <- unit.width * max(reg.np)
	ipl <- .2 # inches per line. Obtained from par->'mai', 'mar'.
	m.bot <- 2; m.lef <- 1; m.top <- 2; m.rig <- 1 # margin size in lines.
	# Convert #gene to image height.
	reg.hei <- sapply(reg.ng, function(r) {
					r * unit.width / pts / reduce.ratio + 
					m.bot * ipl + m.top * ipl  # margins are like intercepts.
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

	plotheat(reg.list, uniq.reg, enrichList, go.algo, ctg.tbl$title, xticks, 
			flood.frac)
	dev.off()
}

# Save plotting data of reads density to a text file.
out.txt <- paste(basename, '.txt', sep='')
out.header <- c(
	'#Do NOT change the following lines if you want to re-draw the image with \
	 replot.r! If you change the matrix values, cut and paste the commented \
	 lines into your new file and run replot.r.',
	paste('#reg2plot:', reg2plot, sep=''),
	paste('#flanksize:', flanksize, sep=''),
	# paste('#intsize:', intsize, sep=''),
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
