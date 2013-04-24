#!/usr/bin/env Rscript
#
# Program: ngs.plot.r
# Purpose: Plot sequencing coverages at different genomic regions.
#          Allow overlaying various coverages with gene lists.
#
# -- by Li Shen, MSSM
# Created:      Nov 2011.
# Last Updated: Mar 2013.
#


# Deal with command line arguments.
cmd.help <- function(){
    cat("\nVisit http://code.google.com/p/ngsplot/wiki/ProgramArguments101 for details\n")
    cat("\nUsage: time.ngs.plot.r -G genome -R region -C [cov|config]file\n")
    cat("                  [Options]\n")
    cat("\n## Mandatory parameters:\n")
    cat("  -G   Genome name, currently supported: mm9, hg19, rn4, sacCer3(ensembl only)\n")
    cat("  -R   Genomic regions to plot: tss, tes, genebody, exon, cgi or *.bed\n")
    cat("  -C   Indexed bam file or a configuration file for multiplot\n")
    cat("  -O   Name for output: multiple files will be generated\n")
    cat("## Optional parameters related to configuration file:\n")
    cat("  -E   Gene list to subset regions\n")
    cat("  -T   Image title\n")
    cat("## Important optional parameters:\n")
    cat("  -F   Further information can be provided to subset regions:\n")
    cat("         for genebody: chipseq(default), rnaseq.\n")
    cat("         for exon: canonical(default), variant, promoter, polyA,\n")
    cat("           altAcceptor, altDonor, altBoth.\n")
    cat("         for cgi: ProximalPromoter(default), Promoter1k, Promoter3k,\n")
    cat("           Genebody, Genedesert, OtherIntergenic, Pericentromere.\n")
    cat("  -D   Gene database: refseq(default), ensembl\n")
    cat("  -I   Shall interval be larger than flanking in plot?(0 or 1, default=automatic)\n")
    cat("  -L   Flanking region size\n")
    cat("  -N   Flanking region factor(will override flanking size)\n")
    cat("  -S   Randomly sample the regions for plot, must be:(0, 1]\n")
    cat("  -P   #CPUs to use. Set 0(default) for auto detection\n")
    # cat("  -PT  #data points to be used in plot(default=100)\n")
    cat("## Misc. parameters:\n")
    cat("  -GO  Gene order algorithm used in heatmaps: total(default), hc, max,\n")
    cat("         prod, diff, pca and none(according to gene list supplied)\n")
    cat("  -CS  Chunk size for loading genes in batch(default=100)\n")
    cat("  -FL  Fragment length used to calculate physical coverage(default=150)\n")
    cat("  -MQ  Mapping quality cutoff to filter reads(default=20)\n")
    cat("  -SE  Shall standard errors be plotted?(0 or 1)\n")
    cat("  -RB  The fraction of extreme values to be trimmed on both ends\n")
    cat("         default=0, 0.05 means 5% of extreme values will be trimmed\n")
    cat("  -FC  Flooding fraction:[0, 1), default=0.02\n")
    cat("  -FI  Forbid image output if set to 1(default=0)\n")
    cat("  -MW  Moving window width to smooth avg. profiles, must be integer\n")
    cat("         1=no(default); 3=slightly; 5=somewhat; 9=quite; 13=super.\n")
    cat("  -H   Opacity of shaded area, suggested value:[0, 0.5]\n")
    cat("         default=0, i.e. no shading, just curves\n")
    # cat("  -E   Calculate weighted coverage according to region size(experimental!)\n")
    cat("\n")
}


###########################################################################
#################### Deal with program input arguments ####################
args <- commandArgs(T)
# args <- unlist(strsplit('-G mm9 -R tss -C H2bub_MB.sorted.bam -O test_full', ' '))

# Program environment variable.
progpath <- Sys.getenv('NGSPLOT')
if(progpath == "") {
    stop("Set environment variable NGSPLOT before run the program. \
          See README for details.\n")
}

# Input argument parser.
source(file.path(progpath, 'lib', 'parse.args.r'))
args.tbl <- parseArgs(args, c('-G', '-C', '-R'))
if(is.null(args.tbl)){
    cmd.help()
    stop('Error in parsing command line arguments. Stop.\n')
}

# Load required libraries.
if(!suppressMessages(require(ShortRead, warn.conflicts=F))) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(ShortRead)
    if(!suppressMessages(require(ShortRead, warn.conflicts=F))) {
        stop('Loading package ShortRead failed!')
    }
}
if(!suppressMessages(require(BSgenome, warn.conflicts=F))) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(BSgenome)
    if(!suppressMessages(require(BSgenome, warn.conflicts=F))) {
        stop('Loading package BSgenome failed!')
    }
}
if(!suppressMessages(require(doMC, warn.conflicts=F))) {
    install.packages('doMC')
    if(!suppressMessages(require(doMC, warn.conflicts=F))) {
        stop('Loading package doMC failed!')
    }
}
if(!suppressMessages(require(caTools, warn.conflicts=F))) {
    install.packages('caTools')
    if(!suppressMessages(require(caTools, warn.conflicts=F))) {
        stop('Loading package caTools failed!')
    }
}
if(!suppressMessages(require(utils, warn.conflicts=F))) {
    install.packages('utils')
    if(!suppressMessages(require(utils, warn.conflicts=F))) {
        stop('Loading package utils failed!')
    }
}

# Configuration: coverage-genelist-title relationships.
ctg.tbl <- ConfigTbl(args.tbl)

# Setup variables from arguments.
argvar.list <- setupVars(args.tbl, ctg.tbl)
genome <- argvar.list$genome  # genome name, such as mm9, hg19, rn4.
reg2plot <- argvar.list$reg2plot  # tss, tes, genebody, bed...
bed.file <- argvar.list$bed.file  # BED file name if reg2plot=bed.
oname <- argvar.list$oname  # output file root name.
fi_tag <- argvar.list$fi_tag  # tag for forbidding image output
lgint <- argvar.list$lgint  # lgint: boolean tag for large interval
flanksize <- argvar.list$flanksize  # flanking region size
flankfactor <- argvar.list$flankfactor  # flanking region factor
samprate <- argvar.list$samprate  # sampling rate
shade.alp <- argvar.list$shade.alp  # shade area alpha
mw <- argvar.list$mw  # moving window width for smooth function.
cores.number <- argvar.list$cores.number  # #CPUs
se <- argvar.list$se  # se: tag for plotting stand errors
robust <- argvar.list$robust  # robust stat fraction
flood.frac <- argvar.list$flood.frac  # flooding fraction.
go.algo <- argvar.list$go.algo  # gene order algorithm used in heatmaps.
gcs <- argvar.list$gcs  # gcs: chunk size for grouping genes.
fraglen <- argvar.list$fraglen  # fragment length for physical coverage.
map.qual <- argvar.list$map.qual  # mapping quality cutoff.
bufsize <- argvar.list$bufsize  # buffer is to ensure smooth coverage at both 
                                # ends of the coverage vector.


# Register doMC with CPU number.
if(cores.number == 0){
    registerDoMC()
} else {
    registerDoMC(cores.number)
}

# Setup plot-related coordinates and variables.
source(file.path(progpath, 'lib', 'genedb.r'))
plotvar.list <- SetupPlotCoord(args.tbl, ctg.tbl, progpath, genome, reg2plot, 
                               lgint, flanksize, bed.file, samprate)
coord.list <- plotvar.list$coord.list  # list of coordinates for unique regions.
rnaseq.gb <- plotvar.list$rnaseq.gb  # tag for RNA-seq data.
lgint <- plotvar.list$lgint  # lgint: automatically determined if not specified.
reg.list <- plotvar.list$reg.list  # region list as in config file.
uniq.reg <- plotvar.list$uniq.reg  # unique region list.
pint <- plotvar.list$pint  # tag for point interval.
exonmodel <- plotvar.list$exonmodel  # exon ranges if rnaseq.gb=True.

# Setup data points for plot.
source(file.path(progpath, 'lib', 'plotlib.r'))
pts.list <- SetPtsSpline(pint, lgint)
pts <- pts.list$pts  # data points for avg. profile and standard errors.
m.pts <- pts.list$m.pts  # middle data points. For pint, m.pts=1.
f.pts <- pts.list$f.pts  # flanking region data points.

# Setup matrix for avg. profiles.
regcovMat <- CreatePlotMat(pts, ctg.tbl)
# Setup matrix for standard errors.
confiMat <- CreateConfiMat(se, pts, ctg.tbl)

# Genomic enrichment for all profiles in the config. Use this for heatmaps.
enrichList <- vector('list', nrow(ctg.tbl))

#######################################################################
# Here start to extract coverages for all genomic regions and calculate 
# data for plotting.

# Load coverage extraction lib.
source(file.path(progpath, 'lib', 'coverage.r'))

# Extract bam file names from configuration and determine if bam-pair is used.
bfl.res <- bamFileList(ctg.tbl)
bam.pair <- bfl.res$bbp  # boolean for bam-pair.
bam.list <- bfl.res$bam.list  # bam file list.

# Determine if bowtie is used to generate the bam files.
# Index bam files if not done yet.
v.map.bowtie <- headerIndexBam(bam.list)

# Retrieve chromosome names for each bam file.
sn.list <- seqnamesBam(bam.list)

# Calculate library size from bam files for normalization.
v.lib.size <- libSizeBam(bam.list)


# ptm <- proc.time()
# Process the config file row by row.
res.time <- system.time(for(r in 1:nrow(ctg.tbl)) {

    reg <- ctg.tbl$glist[r]  # retrieve gene list names.
    # Create coordinate chunk indices.
    chkidx.list <- chunkIndex(nrow(coord.list[[reg]]), gcs)

    # Do coverage for each bam file.
    bam.files <- unlist(strsplit(ctg.tbl$cov[r], ':'))
    # Obtain bam file basic info.
    libsize <- v.lib.size[bam.files[1]]
    sn.inbam <- sn.list[[bam.files[1]]]
    chr.tag <- chrTag(sn.inbam)
    is.bowtie <- v.map.bowtie[bam.files[1]]
    if(class(chr.tag) == 'character') {
        stop(sprintf("Read %s error: %s", bam.files[1], chr.tag))
    }
    result.matrix <- covMatrix(bam.files[1], libsize, sn.inbam, chr.tag, 
                               coord.list[[reg]], chkidx.list, rnaseq.gb, 
                               exonmodel, reg2plot, pint, flanksize, 
                               flankfactor, bufsize, fraglen, map.qual, m.pts, 
                               f.pts, is.bowtie, spit.dot=F)
    if(bam.pair) {  # calculate background.
        pseudo.rpm <- 1e-9
        libsize <- v.lib.size[bam.files[2]]
        sn.inbam <- sn.list[[bam.files[2]]]
        chr.tag <- chrTag(sn.inbam)
        is.bowtie <- v.map.bowtie[bam.files[2]]
        if(class(chr.tag) == 'character') {
            stop(sprintf("Read %s error: %s", bam.files[2], chr.tag))
        }
        bkg.matrix <- covMatrix(bam.files[2], libsize, sn.inbam, chr.tag, 
                                coord.list[[reg]], chkidx.list, rnaseq.gb, 
                                exonmodel, reg2plot, pint, flanksize, 
                                flankfactor, bufsize, fraglen, map.qual, m.pts, 
                                f.pts, is.bowtie, spit.dot=F)
        # browser()
        result.matrix <- log2((result.matrix + pseudo.rpm) / 
                              (bkg.matrix + pseudo.rpm))
    }

    # Calculate SEM if needed. Shut off SEM in single gene case.
    if(nrow(result.matrix) > 1 && se){
        confiMat[, r] <- apply(result.matrix, 2, function(x) CalcSem(x, robust))
    }

    # Book-keep this matrix for heatmap.
    enrichList[[r]] <- result.matrix

    # Return avg. profile.
    regcovMat[, r] <- apply(result.matrix, 2, function(x) mean(x, trim=robust, 
                                                               na.rm=T))
})

print(res.time)




