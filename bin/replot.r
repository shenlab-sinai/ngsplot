#!/usr/bin/env Rscript
#
# Program: replot.r
# Purpose: Re-generate figures using plotting data.
#
# -- by Li Shen, MSSM
# Created:      Nov 2011.
#

# Read environmental variables.
progpath <- Sys.getenv('NGSPLOT')
if(progpath == ""){
    stop("Set environment variable NGSPLOT before run the program. See README for details.\n")
}

source(file.path(progpath, 'lib', 'parse.args.r'))
source(file.path(progpath, 'lib', 'plotlib.r'))

cmd.help <- function(){
    cat("\nVisit http://code.google.com/p/ngsplot/wiki/ProgramArguments101 for details\n")
    cat("\nUsage: replot.r command -I input.zip -O name\n")
    cat("\n  command: prof OR heatmap\n\n")
    cat("## Mandatory parameters:\n")
    cat("  -I   Result zip file created by ngs.plot\n")
    cat("  -O   Output name\n")
    EchoPlotArgs()
    cat("\n")
}

# Read command arguments.
args <- commandArgs(T)
# args <- unlist(strsplit('prof -I test2.zip -O test2.rep -BOX 0 -LEG 0 -VLN 0', ' '))
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

# Parse arguments.
args.tbl <- parseArgs(args, c('-I', '-O'))
if(is.null(args.tbl)){
    cmd.help()
    stop('Error in parsing command line arguments. Stop.\n')
}
iname <- args.tbl['-I']  # iname: input zip file name
oname <- args.tbl['-O']  # oname: output file root name

# Load plotting parameters and data; update plotting parameters.
ifolder <- sub('.zip$', '', basename(iname))
if(ifolder == basename(iname)) {
    stop("Input filename must end with .zip\n")
}
if(command == 'prof') {
    load(unz(iname, file.path(ifolder, 'avgprof.RData')))
    updated.vl <- PlotVars(args.tbl, ls(), prof.misc)
} else if(command == 'heatmap') {
    load(unz(iname, file.path(ifolder, 'heatmap.RData')))
    updated.vl <- PlotVars(args.tbl, ls(), low.count=low.count, 
                           go.paras=go.paras)
    CheckHMColorConfig(hm.color, bam.pair)
} else {
    # pass.
}

################################################
# With updated vars, load plotting procedures and do it.
with(updated.vl, {
    if(command == 'prof') {
        # Update some parameters.
        # Do plot.
        out.plot <- paste(oname, '.pdf', sep='')
        pdf(out.plot, width=plot.width, height=plot.height, 
            pointsize=font.size)
        if(!se) {
            confiMat <- NULL
        }
        if(!suppressMessages(require(caTools, warn.conflicts=F))) {
            install.packages('caTools')
            if(!suppressMessages(require(caTools, warn.conflicts=F))) {
                stop('Loading package caTools failed!')
            }
        }
        plotmat(regcovMat, ctg.tbl$title, ctg.tbl$color, bam.pair, xticks, 
                pts, m.pts, f.pts, pint, shade.alp, confiMat, mw, prof.misc)
        out.dev <- dev.off()
        # Save replot data.
        out.dat <- paste(oname, '.RData', sep='')
        save(plot.width, plot.height, shade.alp, mw, se, prof.misc, 
             file=out.dat)
    } else if(command == 'heatmap') {
        # # Setup multi-core doMC.
        # if(!suppressMessages(require(doMC, warn.conflicts=F))) {
        #     install.packages('doMC')
        #     if(!suppressMessages(require(doMC, warn.conflicts=F))) {
        #         stop('Loading package doMC failed!')
        #     }
        # }
        # # Register doMC with CPU number.
        # registerDoMC(1)

        # Setup heatmap device.
        hd <- SetupHeatmapDevice(reg.list, uniq.reg, ng.list, pts, font.size,
                                 reduce.ratio=rr)
        reg.hei <- hd$reg.hei  # list of image heights for unique regions.
        hm.width <- hd$hm.width  # hm.width: image width.
        hm.height <- hd$hm.height  # hm.height: image height.
        lay.mat <- hd$lay.mat  # lay.mat: matrix for heatmap layout.
        heatmap.mar <- hd$heatmap.mar  # heatmap.mar: heatmap margins.
        out.hm <- paste(oname, '.pdf', sep='')
        pdf(out.hm, width=hm.width, height=hm.height, pointsize=font.size)
        par(mai=heatmap.mar)
        layout(lay.mat, heights=reg.hei)
        # Do heatmap plotting.
        v.low.cutoff <- low.count.ratio * v.low.cutoff
        go.list <- plotheat(reg.list, uniq.reg, enrichList, v.low.cutoff, go.algo, 
                            go.paras, ctg.tbl$title, bam.pair, xticks, flood.frac, 
                            do.plot=T, hm.color=hm.color, color.distr=color.distr, 
                            color.scale=color.scale)
        out.dev <- dev.off()
        # Save replot data.
        out.dat <- paste(oname, '.RData', sep='')
        save(v.low.cutoff, flood.frac, hm.color, go.algo, rr, go.list, 
             color.scale, go.paras, low.count, color.distr, 
             file=out.dat)
    } else {
        # Pass.
    }
})

cat("Done.\n")