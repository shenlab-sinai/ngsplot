#!/usr/bin/env Rscript
#
# Program: ngs.plot.r
# Purpose: Plot sequencing coverages at different genomic regions.
#          Allow overlaying various coverages with gene lists.
#
# -- by Li Shen, MSSM
# Created:      Nov 2011.
# Last Updated: Jul 2013.
#

# library(compiler)
# enableJIT(3)

ngsplot.version <- '2.08'

# Deal with command line arguments.
cmd.help <- function(){
    cat("\nVisit http://code.google.com/p/ngsplot/wiki/ProgramArguments101 for details\n")
    cat(paste("Version:", ngsplot.version, sep=" "))
    cat("\nUsage: ngs.plot.r -G genome -R region -C [cov|config]file\n")
    cat("                  -O name [Options]\n")
    cat("\n## Mandatory parameters:\n")
    cat("  -G   Genome name. Use ngsplotdb.py list to show available genomes.\n")
    cat("  -R   Genomic regions to plot: tss, tes, genebody, exon, cgi, enhancer, dhs or bed\n")
    cat("  -C   Indexed bam file or a configuration file for multiplot\n")
    cat("  -O   Name for output: multiple files will be generated\n")
    cat("## Optional parameters related to configuration file:\n")
    cat("  -E   Gene list to subset regions OR bed file for custom region\n")
    cat("  -T   Image title\n")
    cat("## Important optional parameters:\n")
    cat("  -F   Further information provided to select database table or plottype:\n")
    cat("         This is a string of description separated by comma.\n")
    cat("         E.g. protein_coding,K562,rnaseq(order of descriptors does not matter)\n")
    cat("              means coding genes in K562 cell line drawn in rnaseq mode.\n")
    cat("  -D   Gene database: ensembl(default), refseq\n")
    cat("  -I   Shall interval be larger than flanking in plot?(0 or 1, default=automatic)\n")
    cat("  -L   Flanking region size\n")
    cat("  -N   Flanking region factor(will override flanking size)\n")
    cat("  -S   Randomly sample the regions for plot, must be:(0, 1]\n")
    cat("  -P   #CPUs to use. Set 0(default) for auto detection\n")
    # cat("  -PT  #data points to be used in plot(default=100)\n")
    cat("## Misc. parameters:\n")
    cat("  -GO  Gene order algorithm used in heatmaps: total(default), hc, max,\n")
    cat("         prod, diff, pca and none(according to gene list supplied)\n")
    cat("  -AL  Algorithm used to normalize coverage vectors: spline(default), bin\n")
    cat("  -CS  Chunk size for loading genes in batch(default=100)\n")
    cat("  -FL  Fragment length used to calculate physical coverage(default=150)\n")
    cat("  -MQ  Mapping quality cutoff to filter reads(default=20)\n")
    cat("  -SE  Shall standard errors be plotted?(0 or 1)\n")
    cat("  -RB  The fraction of extreme values to be trimmed on both ends\n")
    cat("         default=0, 0.05 means 5% of extreme values will be trimmed\n")
    cat("  -RZ  Remove all zero profiles in heatmaps(default=1). Set 0 to keep them.\n")
    cat("  -SC  Color scale used to map values to colors in a heatmap.\n")
    cat("         local(default): base on each individual heatmap\n")
    cat("         region: base on all heatmaps belong to the same region\n")
    cat("         global: base on all heatmaps together\n")
    cat("         min_val,max_val: custom scale using a pair of numerics\n")
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
# args <- unlist(strsplit('-G mm9 -R tss -C h3k4me3r.C1.mono.bam -O test', ' '))

# Program environment variable.
progpath <- Sys.getenv('NGSPLOT')
if(progpath == "") {
    stop("Set environment variable NGSPLOT before run the program. \
          See README for details.\n")
}

# Input argument parser.
source(file.path(progpath, 'lib', 'parse.args.r'))
args.tbl <- parseArgs(args, c('-G', '-C', '-R', '-O'))
if(is.null(args.tbl)){
    cmd.help()
    stop('Error in parsing command line arguments. Stop.\n')
}

# Load required libraries.
cat("Loading R libraries")
if(!suppressMessages(require(ShortRead, warn.conflicts=F))) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(ShortRead)
    if(!suppressMessages(require(ShortRead, warn.conflicts=F))) {
        stop('Loading package ShortRead failed!')
    }
}
cat('.')
if(!suppressMessages(require(BSgenome, warn.conflicts=F))) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(BSgenome)
    if(!suppressMessages(require(BSgenome, warn.conflicts=F))) {
        stop('Loading package BSgenome failed!')
    }
}
cat('.')
if(!suppressMessages(require(doMC, warn.conflicts=F))) {
    install.packages('doMC')
    if(!suppressMessages(require(doMC, warn.conflicts=F))) {
        stop('Loading package doMC failed!')
    }
}
cat('.')
if(!suppressMessages(require(caTools, warn.conflicts=F))) {
    install.packages('caTools')
    if(!suppressMessages(require(caTools, warn.conflicts=F))) {
        stop('Loading package caTools failed!')
    }
}
cat('.')
if(!suppressMessages(require(utils, warn.conflicts=F))) {
    install.packages('utils')
    if(!suppressMessages(require(utils, warn.conflicts=F))) {
        stop('Loading package utils failed!')
    }
}
cat('.')
cat("Done\n")

cat("Configuring variables...")
# Load tables of database: default.tbl, dbfile.tbl
default.tbl <- read.delim(file.path(progpath, 'database', 'default.tbl'))
dbfile.tbl <- read.delim(file.path(progpath, 'database', 'dbfile.tbl'))

# Setup variables from arguments.
argvar.list <- setupVars(args.tbl, default.tbl)
genome <- argvar.list$genome  # genome name, such as mm9, hg19, rn4.
reg2plot <- argvar.list$reg2plot  # tss, tes, genebody, bed...
oname <- argvar.list$oname  # output file root name.
galaxy <- argvar.list$galaxy  # tag for Galaxy use.
if(galaxy) {
    avgname <- argvar.list$avgname
    heatmapname <- argvar.list$heatmapname
}
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
color.scale <- argvar.list$color.scale  # string for color scale.
flood.frac <- argvar.list$flood.frac  # flooding fraction.
rm.zero <- argvar.list$rm.zero  # remove all zero tag.
go.algo <- argvar.list$go.algo  # gene order algorithm used in heatmaps.
cov.algo <- argvar.list$cov.algo  # coverage normalization algorithm.
gcs <- argvar.list$gcs  # gcs: chunk size for grouping genes.
fraglen <- argvar.list$fraglen  # fragment length for physical coverage.
map.qual <- argvar.list$map.qual  # mapping quality cutoff.
bufsize <- argvar.list$bufsize  # buffer is to ensure smooth coverage at both 
                                # ends of the coverage vector.

# Configuration: coverage-genelist-title relationships.
ctg.tbl <- ConfigTbl(args.tbl, fraglen)

# Register doMC with CPU number.
if(cores.number == 0){
    registerDoMC()
} else {
    registerDoMC(cores.number)
}

# Setup plot-related coordinates and variables.
source(file.path(progpath, 'lib', 'genedb.r'))
plotvar.list <- SetupPlotCoord(args.tbl, ctg.tbl, default.tbl, dbfile.tbl, 
                               progpath, genome, reg2plot, lgint, flanksize, 
                               samprate, galaxy)
coord.list <- plotvar.list$coord.list  # list of coordinates for unique regions.
rnaseq.gb <- plotvar.list$rnaseq.gb  # tag for RNA-seq data.
lgint <- plotvar.list$lgint  # lgint: automatically determined if not specified.
reg.list <- plotvar.list$reg.list  # region list as in config file.
uniq.reg <- plotvar.list$uniq.reg  # unique region list.
pint <- plotvar.list$pint  # tag for point interval.
exonmodel <- plotvar.list$exonmodel  # exon ranges if rnaseq.gb=True.
Labs <- plotvar.list$Labs # labels for the plot.

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
cat("Done\n")

#######################################################################
# Here start to extract coverages for all genomic regions and calculate 
# data for plotting.

# Load coverage extraction lib.
cat("Analyze bam files and calculate coverage")
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

# Process the config file row by row.
for(r in 1:nrow(ctg.tbl)) {

    reg <- ctg.tbl$glist[r]  # retrieve gene list names.
    # Create coordinate chunk indices.
    chkidx.list <- chunkIndex(nrow(coord.list[[reg]]), gcs)

    # Do coverage for each bam file.
    bam.files <- unlist(strsplit(ctg.tbl$cov[r], ':'))

    # Obtain fraglen for each bam file.
    fraglens <- as.integer(unlist(strsplit(ctg.tbl$fraglen[r], ':')))

    # Obtain bam file basic info.
    libsize <- v.lib.size[bam.files[1]]
    result.pseudo.rpm <- 1e6 / libsize
    sn.inbam <- sn.list[[bam.files[1]]]
    chr.tag <- chrTag(sn.inbam)
    is.bowtie <- v.map.bowtie[bam.files[1]]
    if(class(chr.tag) == 'character') {
        stop(sprintf("Read %s error: %s", bam.files[1], chr.tag))
    }
    # browser()
    result.matrix <- covMatrix(chkidx.list, coord.list[[reg]], rnaseq.gb, 
                               exonmodel, libsize, TRUE, chr.tag, pint, 
                               reg2plot, flanksize, flankfactor, m.pts, f.pts, 
                               bufsize, cov.algo, bam.files[1], sn.inbam, 
                               fraglens[1], map.qual, is.bowtie)
    if(bam.pair) {  # calculate background.
        fraglen2 <- ifelse(length(fraglens) > 1, fraglens[2], fraglens[1])
        libsize <- v.lib.size[bam.files[2]]
        bkg.pseudo.rpm <- 1e6 / libsize
        sn.inbam <- sn.list[[bam.files[2]]]
        chr.tag <- chrTag(sn.inbam)
        is.bowtie <- v.map.bowtie[bam.files[2]]
        if(class(chr.tag) == 'character') {
            stop(sprintf("Read %s error: %s", bam.files[2], chr.tag))
        }
        bkg.matrix <- covMatrix(chkidx.list, coord.list[[reg]], rnaseq.gb, 
                                exonmodel, libsize, TRUE, chr.tag, pint, 
                                reg2plot, flanksize, flankfactor, m.pts, f.pts, 
                                bufsize, cov.algo, bam.files[2], sn.inbam, 
                                fraglen2, map.qual, is.bowtie)
        # browser()
        result.matrix <- log2((result.matrix + result.pseudo.rpm) / 
                              (bkg.matrix + bkg.pseudo.rpm))
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
}
# browser()
cat("Done\n")

########################################
# Add row names to heatmap data.
for(i in 1:length(enrichList)) {
    reg <- ctg.tbl$glist[i]  # gene list name.
    rownames(enrichList[[i]]) <- paste(coord.list[[reg]]$gname, 
                                       coord.list[[reg]]$tid, sep=':')
}
# Some basic parameters.
default.width <- 8  # in inches.
default.height <- 7
xticks <- genXticks(reg2plot, pint, lgint, pts, flanksize, flankfactor, Labs)
unit.width <- 4
rr <- 30  # reduce ratio.
ng.list <- sapply(enrichList, nrow)  # number of genes per heatmap.

# Create image file and plot data into it.
if(!fi_tag){
    cat("Plotting figures...")
    #### Average profile plot. ####
    if(galaxy) {
        out.plot <- avgname
    } else {
        out.plot <- paste(oname, '.avgprof.pdf', sep='')
    }
    pdf(out.plot, width=default.width, height=default.height)
    plotmat(regcovMat, ctg.tbl$title, ctg.tbl$color, bam.pair, xticks, pts, 
            m.pts, f.pts, pint, shade.alp, confiMat, mw)
    out.dev <- dev.off()

    #### Heatmap. ####
    # Setup output device.
    hd <- SetupHeatmapDevice(reg.list, uniq.reg, ng.list, pts, unit.width, rr)
    reg.hei <- hd$reg.hei  # list of image heights for unique regions.
    hm.width <- hd$hm.width  # image width.
    hm.height <- hd$hm.height # image height.
    lay.mat <- hd$lay.mat  # matrix for heatmap layout.
    heatmap.mar <- hd$heatmap.mar # heatmap margins in inches.

    if(galaxy) {
        out.hm <- heatmapname
    } else {
        out.hm <- paste(oname, '.heatmap.pdf', sep='')
    }
    pdf(out.hm, width=hm.width, height=hm.height)
    par(mai=heatmap.mar)
    layout(lay.mat, heights=reg.hei)

    # Do heatmap plotting.
    go.list <- plotheat(reg.list, uniq.reg, enrichList, go.algo, ctg.tbl$title, 
                        bam.pair, xticks, rm.zero, flood.frac, 
                        color.scale=color.scale)
    out.dev <- dev.off()

    cat("Done\n")
} else {
    go.list <- plotheat(reg.list, uniq.reg, enrichList, go.algo, ctg.tbl$title, 
                        bam.pair, xticks, rm.zero, flood.frac, do.plot=F,
                        color.scale=color.scale)
}

# Save plotting data.
if(galaxy==1) { oname1="data" }
cat("Saving results...")
if(galaxy==1) {
   dir.create(oname1, showWarnings=F)
} else {
   dir.create(oname, showWarnings=F)
}
# Average profiles.
if(galaxy==1){
   out.prof <- file.path(oname1, 'avgprof.txt')
}else{
   out.prof <- file.path(oname, 'avgprof.txt')
}
write.table(regcovMat, file=out.prof, row.names=F, sep="\t", quote=F)

# Standard errors of mean.
if(!is.null(confiMat)){
    if(galaxy==1){
       out.confi <- file.path(oname1, 'sem.txt')
    }else{
       out.confi <- file.path(oname, 'sem.txt')
    }
    write.table(confiMat, file=out.confi, row.names=F, sep="\t", quote=F)
}

# Heatmap density values.
for(i in 1:length(enrichList)) {
    reg <- ctg.tbl$glist[i]  # gene list name.
    if(galaxy==1){
       out.heat <- file.path(oname1, paste('hm', i, '.txt', sep=''))
    }else{
       out.heat <- file.path(oname, paste('hm', i, '.txt', sep=''))
    }
    write.table(cbind(coord.list[[reg]][, c('gid', 'gname', 'tid', 'strand')], 
                      enrichList[[i]]), 
                file=out.heat, row.names=F, sep="\t", quote=F)
}

# Avg. profile R data.
if(galaxy==1){
   prof.dat <- file.path(oname1, 'avgprof.RData')
}else{
   prof.dat <- file.path(oname, 'avgprof.RData')
}
save(default.width, default.height, regcovMat, ctg.tbl, bam.pair, xticks, pts, 
     m.pts, f.pts, pint, shade.alp, confiMat, mw, se, file=prof.dat)

# Heatmap R data.
if(galaxy==1){
     heat.dat <- file.path(oname1, 'heatmap.RData')
}else{
     heat.dat <- file.path(oname, 'heatmap.RData')
}
save(reg.list, uniq.reg, ng.list, pts, enrichList, go.algo, ctg.tbl, bam.pair, 
     xticks, rm.zero, flood.frac, unit.width, rr, go.list, color.scale, 
     file=heat.dat)
cat("Done\n")

# Wrap results up.
cat("Wrapping results up...")
cur.dir <- getwd()
if(galaxy==1){
    out.dir <- dirname(oname1)
    out.zip <- basename(oname1)
}else{
    out.dir <- dirname(oname)
    out.zip <- basename(oname)
}
setwd(out.dir)
if(!zip(paste(out.zip, '.zip', sep=''), out.zip, extras='-q')) {
    if(unlink(oname, recursive=T)) {
        warning(sprintf("Unable to delete intermediate result folder: %s", 
                        oname))
    }
}
setwd(cur.dir)
cat("Done\n")
cat("All done. Cheers!\n")

