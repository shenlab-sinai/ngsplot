#!/usr/bin/env Rscript
#
# Program: plotCorrGram.r
# Purpose: Plot CorrGram of the matrix of the output from ngs.plot.
#
# -- by Ning-Yi SHAO, MSSM
# Created:      Jan 2014.
#

cmd.help <- function(){
    cat("\nVisit http://code.google.com/p/ngsplot for details\n")
    cat("\nUsage: plotCorrGram.r -I ngsplot_output.zip [Options]\n")
    cat("\n## Mandatory parameters:\n")
    cat("  -I   Result zip file created by ngs.plot.\n")
    cat("  -O  Output name\n")
    cat("## Optional parameters:\n")
    cat("  -F Functions to be used in calculating the values fed to calculate the correlation (default=mean).\n")
    cat("     mean(default): base on mean of each row.\n")
    cat("     max: base on max of each row.\n")
    cat("     window: base on the center region of each row.\n")
    cat("  -OP Options for -F functions.\n")
    cat("     [0,0.5): For 'mean', trim value for robust estimation, default is 0.\n")
    cat("     window_left, window_right: For 'window', borders of the window (total length of the row is 1), default is 0.33,0.66.\n")
    cat("\n")
}

if(!suppressMessages(require(psych, warn.conflicts=F))) {
        install.packages('psych')
        if(!suppressMessages(require(psych, warn.conflicts=F))) {
            stop('Loading package psych failed!')
        }
    }

if(!suppressMessages(require(corrgram, warn.conflicts=F))) {
        install.packages('corrgram')
        if(!suppressMessages(require(corrgram, warn.conflicts=F))) {
            stop('Loading package corrgram failed!')
        }
    }


# Read command.
args <- commandArgs(T)

if(length(args) < 1) {
    cmd.help()
    stop("No arguments provided.\n")
}

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

#### Check for the function used ####
if('-F' %in% names(args.tbl)){
    stopifnot(args.tbl['-F'] %in% c('mean', 'max', 'window'))
    fun <- args.tbl['-F']
}else{
    fun <- 'mean'
}

#### Check for the options ####
mean.trim <- 0
window.borders <- c(0.33, 0.66)
if('-OP' %in% names(args.tbl)){
    if(fun=='mean'){
        stopifnot(as.numeric(args.tbl['-OP']) >= 0)
        mean.trim <- as.numeric(args.tbl['-OP'])
    }
    if(fun=='window'){
        window.pair <- unlist(strsplit(args.tbl['-OP'], ","))
        if(length(window.pair) != 2 || is.na(as.numeric(window.pair[1])) ||
            is.na(as.numeric(window.pair[2]))) {
                stop("Window borders must be a pair of numerics separated by ','\n")
        }
        window.borders <- c(as.numeric(window.pair[1]), as.numeric(window.pair[2]))
    }
}

# Load plotting parameters and data.
ifolder <- sub('.zip$', '', basename(iname))
if(ifolder == basename(iname)) {
    stop("Input filename must end with .zip\n")
}
load(unz(iname, file.path(ifolder, 'heatmap.RData')))

# getRowSums <- function(x){	
# 	row.sum.x <- rowSums(x)
# 	return(row.sum.x)
# }

getRowMeans <- function(x, mean.trim=0, na.rm=TRUE){
    row.mean.x <- apply(x, 1, mean, trim=mean.trim, na.rm=na.rm)
    names(row.mean.x) <- dimnames(x)[[1]]
    return(row.mean.x)
}

getRowMax <- function(x){  
    row.max.x <- apply(x, 1, max)
    names(row.max.x) <- dimnames(x)[[1]]
    return(row.max.x)
}

getRowWindow <- function(x, window.borders){
    window.x <- x[, c(as.integer(window.borders[1]*length(x[1, ])):as.integer(window.borders[2]*length(x[1, ])))]
    row.window.x <- apply(window.x, 1, mean, na.rm=TRUE)
    names(row.window.x) <- dimnames(x)[[1]]
    return(row.window.x)
}

# calculate the values for the correlation test
switch(fun,
    mean={
        value.list <- lapply(enrichList, getRowMeans, mean.trim=mean.trim)
    },
    max={
        value.list <- lapply(enrichList, getRowMax)
    },
    window={
        value.list <- lapply(enrichList, getRowWindow, window.borders=window.borders)
    }
    )

# prepare the matrix to do correlation test
value.matrix <- t(do.call(mapply, c(cbind, value.list)))
colnames(value.matrix) <- ctg.tbl$title

# do correlation test and write to csv files
pearson <- corr.test(value.matrix, use="pairwise")
write.csv(pearson$r, file=paste(oname, "_pearson_cor.csv", sep=''))
write.csv(pearson$p, file=paste(oname, "_pearson_pvalue.csv", sep=''))

spearman <- corr.test(value.matrix, use="pairwise", method="spearman")
write.csv(spearman$r, file=paste(oname, "_spearman_rho.csv", sep=''))
write.csv(spearman$p, file=paste(oname, "_spearman_pvalue.csv", sep=''))

# plot corrgram
pdf(file=paste(oname, ".pdf", sep=''))
corrgram(value.matrix, upper.panel=panel.pie, lower.panel=panel.ellipse)
dev.off()
