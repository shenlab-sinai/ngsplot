#!/usr/bin/env Rscript  
#
# Program: cov.calc.r
# Purpose: Calculate the genomic coverage for a NGS sample.
# Arguments: genome name, aligned read file, output coverage file.
#
# -- by Li Shen, MSSM
#    Nov 2011
#
# Developers: Li Shen, Ningyi Shao
# Dept. of Neuroscience, Mount Sinai School of Medicine
#

cmd.help <- function(){
	cat("\nUsage: cov.calc.r -G genome_name -R alignedread_file -X align_format -C cov_file [-F frag_len] [-Q mapping_quality] [-V calc_gig_ratio]\n")
	cat("\n")
	cat("Mandatory parameters:\n")
	cat("  -G     Genome name, currently supported: mm9, hg19, rn4.\n")
	cat("  -R     Short read alignment file name(e.g. *.bam, *_export.txt).\n")
	cat("  -X     Alignment format: bam, export(ELAND/CASAVA output).\n")
	cat("  -C     Coverage output file(binary R objects file).\n")
	cat("\nOptional parameters:\n")
	cat("  -P     Number of CPUs to be used(default=1). Set 0 to use all detected CPUs.\n")
	cat("  -F     Expected fragement length(default=100).\n")
	cat("           You need this parameter because our program calculates physical coverage instead of read coverage.\n")
	cat("           Physical coverage makes it handy when comparing samples with different read lengths.\n")
	cat("           It also provides more accurate coverage information.\n")
	cat("  -Q     Mapping quality requirement(default=20).\n")
	cat("           Phred-scaled score which equals -10*log10(probability of error).\n")
	cat("  -V     Calculate genome vs. intragenic ratio of accumulated coverage(default=0, i.e. Off).\n")
	cat("           This ratio can later be used to adjust coveage in plotting.\n")
	cat("           Can be useful when comparing heterochromatic marks with DNA Input.\n")
	cat("\n")
}

args <- commandArgs(T)
progpath <- Sys.getenv('NGSPLOT')
if(progpath == ""){
	stop("Set environment variable NGSPLOT before run the program. See README for details.\n")
}else{
	if(substr(progpath, nchar(progpath), nchar(progpath)) != '/'){	# add trailing slash.
		progpath <- paste(progpath, '/', sep='')
	}
}
source(paste(progpath, 'lib/parse.args.r', sep=''))
args.tbl <- parse.args(args, c('-G', '-R', '-X', '-C'))
if(is.null(args.tbl)){
	cmd.help()
	stop('Error in parsing command line arguments. Stop.\n')
}
genome <- args.tbl['-G']
readfile <- args.tbl['-R']
readformat <- args.tbl['-X']
# Set alignment file format desc string.
if(readformat == 'bam'){
	typestr <- 'BAM'
}else if(readformat == 'export'){
	typestr <- 'SolexaExport'
}else{
	stop('Unsupported alignment file.')
}
covfile <- args.tbl['-C']

if('-F' %in% names(args.tbl)){	# fragment length.
	stopifnot(as.integer(args.tbl['-F']) > 0)
	fraglen <- args.tbl['-F']
}else{
	fraglen <- 100
}
fraglen <- as.integer(fraglen)

if('-Q' %in% names(args.tbl)){	# map quality.
	stopifnot(as.integer(args.tbl['-Q']) > 0)
	mapqual <- args.tbl['-Q']
}else{
	mapqual <- 20
}
mapqual <- as.integer(mapqual)

if('-P' %in% names(args.tbl)){	# set cores number.
	stopifnot(as.integer(args.tbl['-P']) >= 0)
	cores.number <- as.integer(args.tbl['-P'])
}else{
	cores.number <- as.integer(1)
}

# Load chromosome length info for the specified genome.
load(paste(progpath, 'database/genome.size.RData', sep=''))
if(genome == 'mm9'){
	chrlens <- mm9.chrl
}else if(genome == 'hg19'){
	chrlens <- hg19.chrl
}else if(genome == 'rn4'){
	chrlens <- rn4.chrl
}else{
	stop('Unsupported genome name! Exit.')
}

# Load required libraries.
require(ShortRead)||{source("http://bioconductor.org/biocLite.R");biocLite("ShortRead");require(ShortRead)}
require(BSgenome)||{source("http://bioconductor.org/biocLite.R");biocLite("BSgenome");require(BSgenome)}
require(Rsamtools)||{source("http://bioconductor.org/biocLite.R");biocLite("Rsamtools");require(Rsamtools)}
require(doMC)||{install.packages("doMC", dep=T);require(doMC)}

if(cores.number == 0){
	registerDoMC()
} else {
	registerDoMC(cores.number)
}

# Read the big bam file.
NullGR <- function() {
	NullGRanges <- GRanges(seqnames=character(), ranges=IRanges())
	values(NullGRanges) <- DataFrame(score=integer())
	NullGRanges
}

getReads <- function(crn, bam.info, bamfile, mapqual, fraglen, ...){
	need.bam.what <- c('strand', 'pos', 'qwidth', 'mapq')
	which <- RangesList(IRanges(start=0, end=seqlengths(bam.info)[crn]))
	names(which) <- crn
	param <- ScanBamParam(which=which, what=need.bam.what)
	crn <- sub('.fa$', '', crn)
	nochr.i <- grep('^chr', crn, invert=T)
	crn[nochr.i] <- paste('chr', crn[nochr.i], sep='')
	if (!(crn %in% names(chrlens))){
		return(NullGR())
	}
	reads.list <- scanBam(bamfile, param=param)[[1]]
	if (length(reads.list$pos) == 0){
		return(NullGR())
	}else{
	reads.df <- as.data.frame(reads.list)
	if (mapping.method %in% c("TopHat", "Bowtie")){
		reads.df$mapq[is.na(reads.df$mapq)] <- 254
	}
	reads.df <- reads.df[!is.na(reads.df["pos"]),]
	reads.df <- reads.df[reads.df["mapq"] >= mapqual,]
	reads.df["width"] <- fraglen
	reads <- GRanges(seqnames=crn, ranges=IRanges(start=reads.df$pos, width=reads.df$qwidth), score=reads.df$mapq, strand=Rle(reads.df$strand))
    }
}

getAllReads <- function(crns, bamfile, mapqual, fraglen, ...){
	bam.info <- seqinfo(bamfile)
	all.reads <- foreach(k=1:length(crns)) %dopar% {
		getReads(crns[[k]], bam.info, bamfile, mapqual, fraglen, ...)
	}
	suppressWarnings(gc())
	names(all.reads) <- crns
	all.reads
}

# Read alignment file.
source(paste(progpath, 'lib/sep.filename.r', sep=''))
sep.res <- sep.filename(readfile)
if(readformat == 'export'){
	read.aligned <- readAligned(sep.res$path, sep.res$name, type=typestr)
	# Convert chromosome names to convention.
	crn <- chromosome(read.aligned)
	crn <- sub('.fa$', '', crn)
	nochr.i <- grep('^chr', crn, invert=T)
	crn[nochr.i] <- paste('chr', crn[nochr.i], sep='')
	read.aligned <- AlignedRead(sread=sread(read.aligned), id=id(read.aligned), quality=quality(read.aligned),
 		chromosome=as.factor(crn), position=position(read.aligned), strand=strand(read.aligned), 
 		alignQuality=alignQuality(read.aligned))
	# Filter reads aligned to unconventional chromosomes.
	read.filtered <- read.aligned[chromosome(read.aligned) %in% names(chrlens)]
	# Filter reads based on mapping quality.
	read.filtered <- read.filtered[which(quality(alignQuality(read.filtered)) >= mapqual)]	
}else{
	# Read big bam file.
	if(!file.exists(paste(readfile, ".bai", sep=""))){
		indexBam(readfile)
	}
	header <- scanBamHeader(readfile)
	try.mapping <- try(strsplit(header[[1]]$text$'@PG'[[1]], ':')[[1]][2])
	if (class(try.mapping) != "try-error"){
		mapping.method <<- try.mapping
	} else {
		mapping.method <<- "default"
	}
	readfile.bamfile <- BamFile(readfile)
	readfile.info <- seqinfo(readfile.bamfile)
	crns <- seqnames(readfile.info)
	read.aligned.list <- getAllReads(crns, readfile.bamfile, mapqual,
                                 	fraglen)
	read.filtered <- suppressWarnings(do.call(c,
                                          	unname(read.aligned.list)))
	rm(read.aligned.list)
	suppressWarnings(gc())
}

# Calculate genomic coverage.
read.coverage <- coverage(read.filtered, width=chrlens, extend = fraglen - width(read.filtered))
nreads <- length(read.filtered)
mc.reads.coverage <- foreach(k=1:length(read.coverage)) %dopar% {
	read.coverage[[k]]/nreads*1e6
}
names(mc.reads.coverage) <- names(read.coverage)
read.coverage.n <- GenomeData(mc.reads.coverage)
rm(mc.reads.coverage)
suppressWarnings(gc())

# Calculate genome vs. intragenic ratio of coverage.
if('-V' %in% names(args.tbl)){	# fragment length.
	stopifnot(as.integer(args.tbl['-V']) >= 0)
	calc.ratio <- as.integer(args.tbl['-V'])
}else{
	calc.ratio <- 0
}
if(calc.ratio){
	if(genome == 'mm9'){
		load(paste(progpath, 'database/mm9.RData', sep=''))
	}else if(genome == 'hg19'){
		load(paste(progpath, 'database/hg19.RData', sep=''))
	}else if(genome == 'rn4'){
		load(paste(progpath, 'database/rn4.RData', sep=''))
	}else{
		stop('Unsupported genome name! Exit.')
	}
	sum_intragenic_cov <- function(gb_coord, gcov.n){
		cov.sum <- 0
		for(i in 1:nrow(gb_coord)){
			chrn <- gb_coord$chrom[i]
			gb.sta <- ifelse(gb_coord$start[i] <= 1000, 1, gb_coord$start[i]-1000)
			gb.end <- ifelse(gb_coord$end[i]+1000 > length(gcov.n[[chrn]]), length(gcov.n[[chrn]]), gb_coord$end[i]+1000)
			cov.sum <- cov.sum + sum(seqselect(gcov.n[[chrn]], gb.sta, gb.end))
		}
		cov.sum
	}
	cov.sum.refseq <- sum_intragenic_cov(subset(refseq$genebody, byname.uniq), read.coverage.n)
	cov.sum.ensembl <- sum_intragenic_cov(subset(ensembl$genebody, byname.uniq), read.coverage.n)
	cov.sum.genome <- sum(sapply(read.coverage.n, sum))
	gig.ratio.refseq <- cov.sum.genome / cov.sum.refseq
	gig.ratio.ensembl <- cov.sum.genome / cov.sum.ensembl
	# Save coverage data. 
	save(mapqual, readfile, fraglen, genome, nreads, read.coverage, read.coverage.n, gig.ratio.refseq, gig.ratio.ensembl, file=covfile)
}else{
	# Save coverage data. 
	save(mapqual, readfile, fraglen, genome, nreads, read.coverage, read.coverage.n, file=covfile)
}
