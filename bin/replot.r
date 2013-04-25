#!/usr/bin/env Rscript
#
# Program: replot.r
# Purpose: Re-generate figures using plotting data.
#
# -- by Li Shen, MSSM
# Created:      Nov 2011.
# Last Updated: Mar 2013.
#

cmd.help <- function(){
    cat("\nVisit http://code.google.com/p/ngsplot/wiki/ProgramArguments101 for details\n")
    cat("\nUsage: replot.r command -I input.zip -O name\n")
    cat("\n  command: prof OR heatmap\n")
    cat("\n## Mandatory parameters:\n")
    cat("    -I  Result zip file created by ngs.plot\n")
    cat("    -O  Output name\n")
    cat("## Avg. profiles parameters:\n")
    cat("    -SE  Shall standard errors be plotted?(0 or 1)\n")
    cat("    -MW  Moving window width to smooth avg. profiles, must be integer\n")
    cat("           1=no(default); 3=slightly; 5=somewhat; 9=quite; 13=super.\n")
    cat("    -H   Opacity of shaded area, suggested value:[0, 0.5]\n")
    cat("           default=0, i.e. no shading, just curves\n")
    cat("## Heatmap parameters:\n")
    cat("    -GO  Gene order algorithm used in heatmaps: total(default), hc, max,\n")
    cat("           prod, diff, pca and none(according to gene list supplied)\n")
    cat("    -RR  Reduce ratio(default=30). The parameter controls the heatmap height\n")
    cat("           The smaller the value, the taller the heatmap\n")
    cat("    -RZ  Remove all zero profiles in heatmaps(default=1). Set 0 to keep them.\n")
    cat("    -FC  Flooding fraction:[0, 1), default=0.02\n")
    cat("    -P   #CPUs to use. Set 0(default) for auto detection\n")

    # cat("-C     Subset columns to plot using a string such as 1,2,4-6,9.")
    cat("\n")
}

# Read command.
args <- commandArgs(T)
if(length(args) < 1) {
    cmd.help()
    stop("No arguments provided.\n")
}
command <- args[1]
if(!command %in% c('prof', 'heatmap')) {
    stop("Command must be: prof OR heatmap\n")
    cmd.help()
}
args <- args[-1]

# Read environmental variables.
progpath <- Sys.getenv('NGSPLOT')
if(progpath == ""){
    stop("Set environment variable NGSPLOT before run the program. See README for details.\n")
}

# Parse arguments.
source(file.path(progpath, 'lib', 'parse.args.r'))
args.tbl <- parseArgs(args, c('-I', '-O'))
if(is.null(args.tbl)){
    cmd.help()
    stop('Error in parsing command line arguments. Stop.\n')
}
repvar.list <- replotVars(args.tbl)
oname <- repvar.list$oname
iname <- repvar.list$iname
cores.number <- repvar.list$cores.number
# oname: output file root name
# iname: input zip file name
# shade.alp: shade area alpha
# mw: moving window width for smooth function.
# se: tag for plotting stand errors
# flood.frac: flooding fraction.
# go.algo: gene order algorithm used in heatmaps.
# rr: reduce ratio
# cores.number: number of CPUs to use.


################################################
# Load plotting procedures and do it.
ifolder <- sub('.zip$', '', basename(iname))
if(ifolder == basename(iname)) {
    stop("Input filename must end with .zip\n")
}
source(file.path(progpath, 'lib', 'plotlib.r'))

if(command == 'prof') {
    load(unz(iname, file.path(ifolder, 'avgprof.RData')))
    # Update some parameters.
    if('shade.alp' %in% names(repvar.list)) {
        shade.alp <- repvar.list$shade.alp
    }
    if('mw' %in% names(repvar.list)) {
        mw <- repvar.list$mw
    }
    if('se' %in% names(repvar.list)) {
        se <- repvar.list$se
    }
    # Do plot.
    out.plot <- paste(oname, '.pdf', sep='')
    pdf(out.plot, width=default.width, height=default.height)
    if(!se) {
        confiMat <- NULL
    }
    if(!suppressMessages(require(caTools, warn.conflicts=F))) {
        install.packages('caTools')
        if(!suppressMessages(require(caTools, warn.conflicts=F))) {
            stop('Loading package caTools failed!')
        }
    }
    plotmat(regcovMat, ctg.tbl$title, bam.pair, xticks, pts, m.pts, f.pts, pint,
            shade.alp, confiMat, mw)
    dev.off()
} else if(command == 'heatmap') {
    load(unz(iname, file.path(ifolder, 'heatmap.RData')))
    # Setup multi-core doMC.
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

    # Update some parameters.
    if('flood.frac' %in% names(repvar.list)) {
        flood.frac <- repvar.list$flood.frac
    }
    if('rm.zero' %in% names(repvar.list)) {
        rm.zero <- repvar.list$rm.zero
    } else if(!'rm.zero' %in% ls()) {  # backward compatibility.
        rm.zero <- 1
    }
    if('go.algo' %in% names(repvar.list)) {
        go.algo <- repvar.list$go.algo
    }
    if('rr' %in% names(repvar.list)) {
        rr <- repvar.list$rr
    }
    # Setup heatmap device.
    hd <- SetupHeatmapDevice(reg.list, uniq.reg, ng.list, pts, reduce.ratio=rr)
    reg.hei <- hd$reg.hei  # list of image heights for unique regions.
    hm.width <- hd$hm.width  # hm.width: image width.
    hm.height <- hd$hm.height  # hm.height: image height.
    lay.mat <- hd$lay.mat  # lay.mat: matrix for heatmap layout.
    heatmap.mar <- hd$heatmap.mar  # heatmap.mar: heatmap margins.
    # Do plot.
    out.hm <- paste(oname, '.pdf', sep='')
    pdf(out.hm, width=hm.width, height=hm.height)
    par(mar=heatmap.mar)
    layout(lay.mat, heights=reg.hei)
    # Do heatmap plotting.
    plotheat(reg.list, uniq.reg, enrichList, go.algo, ctg.tbl$title, bam.pair, 
             xticks, rm.zero, flood.frac)
    dev.off()
} else {
    # Pass.
}

cat("Done.\n")