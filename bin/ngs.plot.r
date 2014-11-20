#!/usr/bin/env Rscript
#
# Program: ngs.plot.r
# Purpose: Plot sequencing coverages at different genomic regions.
#          Allow overlaying various coverages with gene lists.
#
# -- by Li Shen, MSSM
# Created:      Nov 2011.
#

# library(compiler)
# enableJIT(3)

ngsplot.version <- '2.47'
# Program environment variable.
progpath <- Sys.getenv('NGSPLOT')
if(progpath == "") {
    stop("Set environment variable NGSPLOT before run the program. See README 
for details.\n")
}

source(file.path(progpath, 'lib', 'parse.args.r'))
source(file.path(progpath, 'lib', 'genedb.r'))
source(file.path(progpath, 'lib', 'plotlib.r'))
source(file.path(progpath, 'lib', 'coverage.r'))

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
    EchoCoverageArgs()
    EchoPlotArgs()
    cat("\n")
}


###########################################################################
#################### Deal with program input arguments ####################
args <- commandArgs(T)
# args <- unlist(strsplit('-G mm9 -R tss -C WT_IN_SORTED.bam -O WT_IN -Debug 1', ' '))

# Input argument parser.
args.tbl <- parseArgs(args, c('-G', '-C', '-R', '-O'))
if(is.null(args.tbl)){
    cmd.help()
    stop('Error in parsing command line arguments. Stop.\n')
}
genome <- args.tbl['-G']
reg2plot <- args.tbl['-R']
oname <- args.tbl['-O']

cat("Configuring variables...")
# Load tables of database: default.tbl, dbfile.tbl
default.tbl <- read.delim(file.path(progpath, 'database', 'default.tbl'))
dbfile.tbl <- read.delim(file.path(progpath, 'database', 'dbfile.tbl'))
CheckRegionAllowed(reg2plot, default.tbl)

# Setup variables from arguments.
cov.args <- CoverageVars(args.tbl, reg2plot)
attach(cov.args)
plot.args <- PlotVars(args.tbl)
attach(plot.args)

# Configuration: coverage-genelist-title relationships.
ctg.tbl <- ConfigTbl(args.tbl, fraglen)

# Setup plot-related coordinates and variables.
plotvar.list <- SetupPlotCoord(args.tbl, ctg.tbl, default.tbl, dbfile.tbl, 
                               progpath, genome, reg2plot, inttag, flanksize, 
                               samprate, galaxy)
attach(plotvar.list)

# Setup data points for plot.
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
# Register doMC with CPU number.
if(cores.number == 0){
    registerDoMC()
} else {
    registerDoMC(cores.number)
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

#######################################################################
# Here start to extract coverages for all genomic regions and calculate 
# data for plotting.

cat("Analyze bam files and calculate coverage")
# Extract bam file names from configuration and determine if bam-pair is used.
bfl.res <- bamFileList(ctg.tbl)
bam.pair <- bfl.res$bbp  # boolean for bam-pair.
bam.list <- bfl.res$bam.list  # bam file list.
CheckHMColorConfig(hm.color, bam.pair)

# Determine if bowtie is used to generate the bam files.
# Index bam files if not done yet.
v.map.bowtie <- headerIndexBam(bam.list)

# Retrieve chromosome names for each bam file.
sn.list <- seqnamesBam(bam.list)

# Calculate library size from bam files for normalization.
v.lib.size <- libSizeBam(bam.list)

v.low.cutoff <- vector("integer", nrow(ctg.tbl))  # low count cutoffs.
# Process the config file row by row.
# browser()
for(r in 1:nrow(ctg.tbl)) {  # r: index of plots/profiles.

    reg <- ctg.tbl$glist[r]  # retrieve gene list names.
    # Create coordinate chunk indices.
    chkidx.list <- chunkIndex(nrow(coord.list[[reg]]), gcs)

    # Do coverage for each bam file.
    bam.files <- unlist(strsplit(ctg.tbl$cov[r], ':'))

    # Obtain fraglen for each bam file.
    fraglens <- as.integer(unlist(strsplit(ctg.tbl$fraglen[r], ':')))

    # Obtain bam file basic info.
    libsize <- v.lib.size[bam.files[1]]
    v.low.cutoff[r] <- low.count / libsize * 1e6
    result.pseudo.rpm <- 1e6 / libsize
    sn.inbam <- sn.list[[bam.files[1]]]
    chr.tag <- chrTag(sn.inbam)
    is.bowtie <- v.map.bowtie[bam.files[1]]
    if(class(chr.tag) == 'character') {
        stop(sprintf("Read %s error: %s", bam.files[1], chr.tag))
    }
    # browser()
    # Rprof("Rprof_covBamExons2.out", append=T)
    result.matrix <- covMatrix(debug, chkidx.list, coord.list[[reg]], rnaseq.gb, 
                               exonmodel, libsize, TRUE, chr.tag, pint, 
                               reg2plot, flanksize, flankfactor, m.pts, f.pts, 
                               bufsize, cov.algo, bam.files[1], sn.inbam, 
                               fraglens[1], map.qual, is.bowtie, 
                               strand.spec=strand.spec)
    # Rprof(NULL)
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
        bkg.matrix <- covMatrix(debug, chkidx.list, coord.list[[reg]], rnaseq.gb, 
                                exonmodel, libsize, TRUE, chr.tag, pint, 
                                reg2plot, flanksize, flankfactor, m.pts, f.pts, 
                                bufsize, cov.algo, bam.files[2], sn.inbam, 
                                fraglen2, map.qual, is.bowtie, 
                                strand.spec=strand.spec)
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
xticks <- genXticks(reg2plot, pint, lgint, pts, flanksize, flankfactor, Labs)
unit.width <- 4
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
    pdf(out.plot, width=plot.width, height=plot.height, pointsize=font.size)
    plotmat(regcovMat, ctg.tbl$title, ctg.tbl$color, bam.pair, xticks, pts, 
            m.pts, f.pts, pint, shade.alp, confiMat, mw, prof.misc)
    out.dev <- dev.off()

    #### Heatmap. ####
    # Setup output device.
    hd <- SetupHeatmapDevice(reg.list, uniq.reg, ng.list, pts, font.size, 
                             unit.width, rr)
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
    pdf(out.hm, width=hm.width, height=hm.height, pointsize=font.size)
    par(mai=heatmap.mar)
    layout(lay.mat, heights=reg.hei)

    # Do heatmap plotting.
    go.list <- plotheat(reg.list, uniq.reg, enrichList, v.low.cutoff, go.algo, 
                        go.paras, ctg.tbl$title, bam.pair, xticks, flood.frac, 
                        do.plot=T, hm.color=hm.color, color.distr=color.distr, 
                        color.scale=color.scale)
    out.dev <- dev.off()
    cat("Done\n")
} else {
    go.list <- plotheat(reg.list, uniq.reg, enrichList, v.low.cutoff, go.algo, 
                        go.paras, ctg.tbl$title, bam.pair, xticks, flood.frac, 
                        do.plot=F, hm.color=hm.color, color.distr=color.distr, 
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
save(plot.width, plot.height, regcovMat, ctg.tbl, bam.pair, xticks, pts, 
     m.pts, f.pts, pint, shade.alp, confiMat, mw, prof.misc, se, v.lib.size, 
     font.size,
     file=prof.dat)

# Heatmap R data.
if(galaxy==1) {
     heat.dat <- file.path(oname1, 'heatmap.RData')
} else {
     heat.dat <- file.path(oname, 'heatmap.RData')
}
save(reg.list, uniq.reg, ng.list, pts, enrichList, v.low.cutoff, go.algo, 
     ctg.tbl, bam.pair, xticks, flood.frac, hm.color, unit.width, rr, 
     go.list, color.scale, v.lib.size, font.size, go.paras, low.count,
     color.distr, 
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

