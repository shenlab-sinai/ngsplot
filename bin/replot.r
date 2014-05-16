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
    cat("## Optional parameters:\n")
    cat("    -FS Font size(default=12).\n")
    cat("## Avg. profiles parameters:\n")
    cat("    -SE  Shall standard errors be plotted?(0 or 1)\n")
    cat("    -MW  Moving window width to smooth avg. profiles, must be integer\n")
    cat("           1=no(default); 3=slightly; 5=somewhat; 9=quite; 13=super.\n")
    cat("    -LEG Draw legend? 1(default) or 0.\n")
    cat("    -BOX Draw box around plot? 1(default) or 0.\n")
    cat("    -VLN Draw vertical lines? 1(default) or 0.\n")
    cat("    -XYL Draw X- and Y-axis labels? 1(default) or 0.\n")
    cat("    -LWD Line width(default=3).\n")
    cat("    -H   Opacity of shaded area, suggested value:[0, 0.5]\n")
    cat("           default=0, i.e. no shading, just curves\n")
    cat("## Heatmap parameters:\n")
    cat("    -GO  Gene order algorithm used in heatmaps: total(default), hc, max,\n")
    cat("           prod, diff, pca and none(according to gene list supplied)\n")
    cat("    -RR  Reduce ratio(default=30). The parameter controls the heatmap height\n")
    cat("           The smaller the value, the taller the heatmap\n")
    cat("    -RZ  Remove all zero profiles in heatmaps(default=1). Set 0 to keep them.\n")
    cat("    -SC  Color scale used to map values to colors in a heatmap.\n")
    cat("           local(default): base on each individual heatmap\n")
    cat("           region: base on all heatmaps belong to the same region\n")
    cat("           global: base on all heatmaps together\n")
    cat("           min_val,max_val: custom scale using a pair of numerics\n")
    cat("    -FC  Flooding fraction:[0, 1), default=0.02\n")
    cat("    -CO  Color for heatmap. For bam-pair, use color-pair(neg_color:pos_color).\n")
    cat("           Hint: must use R colors, such as darkgreen, yellow and blue2.\n")
    cat("    -P   #CPUs to use. Set 0(default) for auto detection\n")

    # cat("-C     Subset columns to plot using a string such as 1,2,4-6,9.")
    cat("\n")
}

# Read command.
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
iname <- args.tbl['-I']  # iname: input zip file name
oname <- args.tbl['-O']  # oname: output file root name

# Load plotting parameters and data.
ifolder <- sub('.zip$', '', basename(iname))
if(ifolder == basename(iname)) {
    stop("Input filename must end with .zip\n")
}
if(command == 'prof') {
    load(unz(iname, file.path(ifolder, 'avgprof.RData')))
} else if(command == 'heatmap') {
    load(unz(iname, file.path(ifolder, 'heatmap.RData')))
} else {
    # pass.
}
# Update or add new variables to the environment.
existing.vl <- ls()
repvar.list <- replotVars(args.tbl, existing.vl, bam.pair, .go.allowed)
for(i in 1:length(repvar.list)) {
    vn <- names(repvar.list)[i]
    eval(parse(text=sprintf("%s <- repvar.list$%s", vn, vn)))
}
if(!"color" %in% colnames(ctg.tbl)) {  # backward compatibility.
    ctg.tbl <- data.frame(ctg.tbl, color=NA)
}


################################################
# Load plotting procedures and do it.
source(file.path(progpath, 'lib', 'plotlib.r'))
if(command == 'prof') {
    # load(unz(iname, file.path(ifolder, 'avgprof.RData')))
    # Update some parameters.
    # Do plot.
    out.plot <- paste(oname, '.pdf', sep='')
    pdf(out.plot, width=default.width, height=default.height, 
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
    save(shade.alp, mw, se, file=out.dat)

} else if(command == 'heatmap') {
    # load(unz(iname, file.path(ifolder, 'heatmap.RData')))
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

    # Setup heatmap device.
    hd <- SetupHeatmapDevice(reg.list, uniq.reg, ng.list, pts, reduce.ratio=rr)
    reg.hei <- hd$reg.hei  # list of image heights for unique regions.
    hm.width <- hd$hm.width  # hm.width: image width.
    hm.height <- hd$hm.height  # hm.height: image height.
    lay.mat <- hd$lay.mat  # lay.mat: matrix for heatmap layout.
    heatmap.mar <- hd$heatmap.mar  # heatmap.mar: heatmap margins.
    # Do plot.
    out.hm <- paste(oname, '.pdf', sep='')
    pdf(out.hm, width=hm.width, height=hm.height, pointsize=font.size)
    par(mai=heatmap.mar)
    layout(lay.mat, heights=reg.hei)
    # Do heatmap plotting.
    go.list <- plotheat(reg.list, uniq.reg, enrichList, v.low.cutoff, go.algo, ctg.tbl$title, 
                        bam.pair, xticks, rm.zero, flood.frac, 
                        hm.color=hm.color, color.scale=color.scale)
    out.dev <- dev.off()
    # Save replot data.
    out.dat <- paste(oname, '.RData', sep='')
    save(flood.frac, rm.zero, go.algo, rr, go.list, file=out.dat)
} else {
    # Pass.
}

cat("Done.\n")