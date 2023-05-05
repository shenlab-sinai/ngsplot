#### plotlib.r ####
# This contains the library for plotting related functions.
#
# Authors: Li Shen, Ningyi Shao
# 
# Created: Feb 19, 2013
# Last updated: May 21, 2013
#

SetupHeatmapDevice <- function(reg.list, uniq.reg, ng.list, pts, font.size=12,
                               unit.width=4, reduce.ratio=30) {
# Configure parameters for heatmap output device. The output is used by 
# external procedures to setup pdf device ready for heatmap plotting.
# Args:
#   reg.list: region list as in config file.
#   uniq.reg: unique region list.
#   ng.list: number of genes per heatmap in the order as config file.
#   pts: data points (number of columns of heatmaps).
#   font.size: font size.
#   unit.width: image width per heatmap.
#   reduce.ratio: how compressed are genes in comparison to data points? This 
#                 controls image height.

    # Number of plots per region.
    reg.np <- sapply(uniq.reg, function(r) sum(reg.list==r))

    # Number of genes per region.
    reg.ng <- sapply(uniq.reg, function(r) {
                    ri <- which(reg.list==r)[1]
                    ng.list[ri]
                })

    # Adjustment ratio.
    origin.fs <- 12  # default font size.
    fs.adj.ratio <- font.size / origin.fs
    # Margin size (in lines) adjusted by ratio.
    m.bot <- fs.adj.ratio * 2
    m.lef <- fs.adj.ratio * 1.5
    m.top <- fs.adj.ratio * 2
    m.rig <- fs.adj.ratio * 1.5 
    key.in <- fs.adj.ratio * 1.0  # colorkey in inches.
    m.lef.diff <- (fs.adj.ratio - 1) * 1.5
    m.rig.diff <- (fs.adj.ratio - 1) * 1.5
    # Setup image size.
    hm.width <- (unit.width + m.lef.diff + m.rig.diff) * max(reg.np)
    ipl <- .2 # inches per line. Obtained from par->'mai', 'mar'.
    # Convert #gene to image height.
    reg.hei <- sapply(reg.ng, function(r) {
                    c(key.in,  # colorkey + margin.
                      r * unit.width / pts / reduce.ratio + 
                      m.bot * ipl + m.top * ipl)  # heatmap + margin.
                })
    reg.hei <- c(reg.hei)
    hm.height <- sum(reg.hei)

    # Setup layout of the heatmaps.
    lay.mat <- matrix(0, ncol=max(reg.np), nrow=length(reg.np) * 2)
    fig.n <- 1  # figure plotting number.
    for(i in 1:length(reg.np)) {
        row.upper <- i * 2 - 1
        row.lower <- i * 2
        for(j in 1:reg.np[i]) {
            lay.mat[row.upper, j] <- fig.n;
            fig.n <- fig.n + 1
            lay.mat[row.lower, j] <- fig.n;
            fig.n <- fig.n + 1
        }
    }

    list(reg.hei=reg.hei, hm.width=hm.width, hm.height=hm.height, 
         lay.mat=lay.mat, heatmap.mar=c(m.bot, m.lef, m.top, m.rig) * ipl)
}

SetPtsSpline <- function(pint, lgint) {
# Set data points for spline function.
# Args:
#   pint: tag for point interval.
# Return: list of data points, middle data points, flanking data points.

    pts <- 100  # data points to plot: 0...pts
    if(pint){  # point interval.
        m.pts <- 1  # middle part points.
        f.pts <- pts / 2  # flanking part points.
    } else {
        if(lgint) {
            m.pts <- pts / 5 * 3 + 1
            f.pts <- pts / 5 + 1
        } else {
            m.pts <- pts / 5 + 1
            f.pts <- pts / 5 * 2 + 1
        }
    }
    list(pts=pts, m.pts=m.pts, f.pts=f.pts)
}

CreatePlotMat <- function(pts, ctg.tbl) {
# Create matrix for avg. profiles.
# Args:
#   pts: data points.
#   ctg.tbl: configuration table.
# Return: avg. profile matrix initialized to zero.

    regcovMat <- matrix(0, nrow=pts + 1, ncol=nrow(ctg.tbl))
    colnames(regcovMat) <- ctg.tbl$title
    regcovMat
}

CreateConfiMat <- function(se, pts, ctg.tbl){
# Create matrix for standard errors.
# Args:
#   se: tag for standard error plotting.
#   pts: data points.
#   ctg.tbl: configuration table.
# Return: standard error matrix initialized to zero or null.

    if(se){
        confiMat <- matrix(0, nrow=pts + 1, ncol=nrow(ctg.tbl))
        colnames(confiMat) <- ctg.tbl$title
    } else {
        confiMat <- NULL
    }
    confiMat
}

col2alpha <- function(col2use, alpha){
# Convert a vector of solid colors to semi-transparent colors.
# Args:
#   col2use: vector of colors.
#   alpha: represents degree of opacity - [0,1]
# Return: vector of transformed colors.

    apply(col2rgb(col2use), 2, function(x){
            rgb(x[1], x[2], x[3], alpha=alpha*255, maxColorValue=255)
        })
}

smoothvec <- function(v, radius, method=c('mean', 'median')){
# Given a vector of coverage, return smoothed version of coverage.
# Args:
#   v: vector of coverage
#   radius: fraction of org. vector size.
#   method: smooth method
# Return: vector of smoothed coverage.

    stopifnot(is.vector(v))
    stopifnot(length(v) > 0)
    stopifnot(radius > 0 && radius < 1)

    halfwin <- ceiling(length(v) * radius)
    s <- rep(NA, length(v))

    for(i in 1:length(v)){
        winpos <- (i - halfwin) : (i + halfwin)
        winpos <- winpos[winpos > 0 & winpos <= length(v)]
        if(method == 'mean'){
            s[i] <- mean(v[winpos])
        }else if(method == 'median'){
            s[i] <- median(v[winpos])
        }
    }
    s
}

smoothplot <- function(m, radius, method=c('mean', 'median')){
# Smooth the entire avg. profile matrix using smoothvec.
# Args:
#   m: avg. profile matrix
#   radius: fraction of org. vector size.
#   method: smooth method.
# Return: smoothed matrix.

    stopifnot(is.matrix(m) || is.vector(m))

    if(is.matrix(m)) {
        for(i in 1:ncol(m)) {
            m[, i] <- smoothvec(m[, i], radius, method)
        }
    } else {
        m <- smoothvec(m, radius, method)
    }
    m
}

genXticks <- function(reg2plot, pint, lgint, pts, flanksize, flankfactor, 
                      Labs) {
# Generate X-ticks for plotting.
# Args:
#   reg2plot: string representation of region.
#   pint: point interval.
#   lgint: tag for large interval.
#   pts: data points.
#   flanksize: flanking region size in bps.
#   flankfactor: flanking region factor.
#   Labs: character vector of labels of the genomic region.
# Return: list of x-tick position and label.

    if(pint){   # point interval.
        mid.lab <- Labs[1]
        tick.pos <- c(0, pts / 4, pts / 2, pts / 4 * 3, pts)
        tick.lab <- as.character(c(-flanksize, -flanksize/2, mid.lab, 
            flanksize/2, flanksize))
    }else{
        left.lab <- Labs[1]
        right.lab <- Labs[2]
        tick.pos <- c(0, pts / 5, pts / 5 * 2, pts / 5 * 3, pts / 5 * 4, pts)
        if(lgint){  # large interval: fla int int int fla
            if(flankfactor > 0){  # show percentage at x-tick.
                tick.lab <- c(sprintf("%d%%", -flankfactor*100), 
                              left.lab, '33%', '66%', right.lab,
                              sprintf("%d%%", flankfactor*100))
            } else{  # show bps at x-tick.
                tick.lab <- c(as.character(-flanksize),
                              left.lab, '33%', '66%', right.lab,
                              as.character(flanksize))
            }
        } else {    # small interval: fla fla int fla fla.
            if(flankfactor > 0){
                tick.lab <- c(sprintf("%d%%", -flankfactor*100), 
                              sprintf("%d%%", -flankfactor*50), 
                              left.lab, right.lab,
                              sprintf("%d%%", flankfactor*50), 
                              sprintf("%d%%", flankfactor*100))
            } else {
                tick.lab <- c(as.character(-flanksize),
                              as.character(-flanksize/2),
                              left.lab, right.lab,
                              as.character(flanksize/2),
                              as.character(flanksize))
            }
        }
    }
    list(pos=tick.pos, lab=tick.lab)
}

plotmat <- function(regcovMat, title2plot, plot.colors, bam.pair, xticks, 
                    pts, m.pts, f.pts, pint, shade.alp=0, confiMat=NULL, mw=1, 
                    misc.options=list(yscale='auto', legend=T, box=T, vline=T, 
                                      xylab=T, line.wd=3)) {
# Plot avg. profiles and standard errors around them.
# Args:
#   regcovMat: matrix for avg. profiles.
#   title2plot: profile names, will be shown in figure legend.
#   plot.colors: vector of color specifications for all curves.
#   bam.pair: boolean for bam-pair data.
#   xticks: X-axis ticks.
#   pts: data points
#   m.pts: middle part data points
#   f.pts: flanking part data points
#   pint: tag for point interval
#   shade.alp: shading area alpha
#   confiMat: matrix for standard errors.
#   mw: moving window size for smoothing function.
#   misc.options: list of misc. options - y-axis scale, legend, box around plot, 
#       verticle lines, X- and Y-axis labels, line width.

    # Smooth avg. profiles if specified.
    if(mw > 1){
        regcovMat <- as.matrix(runmean(regcovMat, k=mw, alg='C', 
                                       endrule='mean'))
    }

    # Choose colors.
    if(any(is.na(plot.colors))) {
        ncurve <- ncol(regcovMat)
        if(ncurve <= 8) {
            suppressMessages(require(RColorBrewer, warn.conflicts=F))
            col2use <- brewer.pal(ifelse(ncurve >= 3, ncurve, 3), 'Dark2')
            col2use <- col2use[1:ncurve]
        } else {
            col2use <- rainbow(ncurve)
        }
    } else {
        col2use <- plot.colors
    }
    col2use <- col2alpha(col2use, 0.8)

    # Plot profiles.
    ytext <- ifelse(bam.pair, "log2(Fold change vs. control)", 
                              "Read count Per Million mapped reads")
    xrange <- 0:pts
    y.lim <- NULL
    if(length(misc.options$yscale) == 2) {
        y.lim <- misc.options$yscale
    }
    matplot(xrange, regcovMat, 
            xaxt='n', type="l", col=col2use, ylim=y.lim,
            lty="solid", lwd=misc.options$line.wd, frame.plot=F, ann=F)
    if(misc.options$xylab) {
        title(xlab="Genomic Region (5' -> 3')", ylab=ytext)
    }
    axis(1, at=xticks$pos, labels=xticks$lab, lwd=3, lwd.ticks=3)
    if(misc.options$box) {
        # box around plot.
        box()
    }

    # Add shade area.
    if(shade.alp > 0){
        for(i in 1:ncol(regcovMat)){
            v.x <- c(xrange[1], xrange, xrange[length(xrange)])
            v.y <- regcovMat[, i]
            v.y <- c(0, v.y, 0)
            col.rgb <- col2rgb(col2use[i])
            p.col <- rgb(col.rgb[1, 1], col.rgb[2, 1], col.rgb[3, 1], 
                         alpha=shade.alp * 255, maxColorValue=255)
            polygon(v.x, v.y, density=-1, border=NA, col=p.col)
        }
    }

    # Add standard errors.
    if(!is.null(confiMat)){
        v.x <- c(xrange, rev(xrange))
        for(i in 1:ncol(confiMat)){
            v.y <- c(regcovMat[, i] + confiMat[, i], 
                     rev(regcovMat[, i] - confiMat[, i]))
            col.rgb <- col2rgb(col2use[i])
            p.col <- rgb(col.rgb[1, 1], col.rgb[2, 1], col.rgb[3, 1], 
                         alpha=0.2 * 255, maxColorValue=255)
            polygon(v.x, v.y, density=-1, border=NA, col=p.col)
        }
    }

    if(misc.options$vline) {
        # Add gray lines indicating feature boundaries.
        if(pint) {
            abline(v=f.pts, col="gray", lwd=2)
        } else {
            abline(v=f.pts - 1, col="gray", lwd=2)
            abline(v=f.pts + m.pts - 2, col="gray", lwd=2)
        }
    }

    if(misc.options$legend) {
        # Legend.
        legend("topright", title2plot, text.col=col2use)
    }
}

spline_mat <- function(mat, n=100){
# Calculate splined coverage for a matrix.
# Args:
#   mat: each column represents a profile to be interpolated.
#   n: number of data points to be interpolated.

    foreach(r=iter(mat, by='row'), 
        .combine='rbind', .multicombine=T) %dopar% {
        spline(1:length(r), r, n)$y
    }
}

RankNormalizeMatrix <- function(mat, low.cutoff) {
# Rank-based normalization for a data matrix.
# Args:
#   mat: data matrix.
#   low.cutoff: low value cutoff.
# Return: rank normalized matrix.

    stopifnot(is.matrix(mat))

    concat.dat <- c(mat)
    low.mask <- concat.dat < low.cutoff
    concat.r <- rank(concat.dat)
    concat.r[low.mask] <- 0

    matrix(concat.r, nrow=nrow(mat))
}

OrderGenesHeatmap <- function(enrichList, lowCutoffs,
                              method=c('total', 'max', 'prod', 'diff', 'hc', 
                                       'none', 'km'),
                              go.paras=list(knc=5, max.iter=20, nrs=30)) {
# Order genes with a list of heatmap data.
# Args: 
#   enrichList: heatmap data in a list.
#   lowCutoffs: low count cutoff for normalized count data.
#   method: algorithm used to order genes.
#   go.paras: gene ordering parameters.
# Returns: a vector of REVERSED gene orders. 
# NOTE: due to the design of image function, the order needs to be reversed
#       so that the genes will be shown correctly. 

    rankList <- mapply(RankNormalizeMatrix, 
                       mat=enrichList, low.cutoff=lowCutoffs, SIMPLIFY=F)
    np <- length(enrichList)

    if(method == 'hc') {
        rankCombined <- do.call('cbind', rankList)
        # Clustering and order genes.
        hc <- hclust(dist(rankCombined, method='euclidean'), 
                     method='complete')
        memb <- cutree(hc, k = go.paras$knc)
        list(rev(hc$order), memb)  # reversed!
    } else if(method == 'km') {
        rankCombined <- do.call('cbind', rankList)
        km <- kmeans(rankCombined, centers=go.paras$knc, 
                     iter.max=go.paras$max.iter, nstart=go.paras$nrs)
        list(rev(order(km$cluster)), km$cluster)  # reversed!
    } else if(method == 'total' || method == 'diff' && np == 1) {  
        list(order(rowSums(rankList[[1]])), NULL)
    } else if(method == 'max') {  # peak enrichment value of the 1st profile.
        list(order(apply(rankList[[1]], 1, max)), NULL)
    } else if(method == 'prod') {  # product of all profiles.
        rs.mat <- sapply(rankList, rowSums)
        g.prod <- apply(rs.mat, 1, prod)
        list(order(g.prod), NULL)
    } else if(method == 'diff' && np > 1) {  # difference between 2 profiles.
        list(order(rowSums(rankList[[1]]) - rowSums(rankList[[2]])), NULL)
    } else if(method == 'none') {  # according to the order of input gene list.
        # Because the image function draws from bottom to top, the rows are 
        # reversed to give a more natural look.
        list(rev(1:nrow(enrichList[[1]])),NULL)
    } else {
        # pass.
    }
}


plotheat <- function(reg.list, uniq.reg, enrichList, v.low.cutoff, go.algo, 
                     go.paras, title2plot, bam.pair, xticks, flood.q=.02, 
                     do.plot=T, hm.color="default", color.distr=.6, 
                     color.scale='local') {
# Plot heatmaps with genes ordered according to some algorithm.
# Args:
#   reg.list: factor vector of regions as in configuration.
#   uniq.reg: character vector of unique regions.
#   enrichList: list of heatmap data.
#   v.low.cutoff: low count cutoff for normalized count data.
#   go.algo: gene order algorithm.
#   go.paras: gene ordering parameters.
#   title2plot: title for each heatmap. Same as the legends in avgprof.
#   bam.pair: boolean tag for bam-pair.
#   xticks: info for X-axis ticks.
#   flood.q: flooding percentage.
#   do.plot: boolean tag for plotting heatmaps.
#   hm.color: string for heatmap colors.
#   color.distr: positive number for color distribution.
#   color.scale: string for the method to adjust color scale.
# Returns: ordered gene names for each unique region as a list.

    # Setup basic parameters.
    ncolor <- 256
    if(bam.pair) {
        if(hm.color != "default") {
            tri.colors <- unlist(strsplit(hm.color, ':'))
            neg.color <- tri.colors[1]
            if(length(tri.colors) == 2) {
                neu.color <- 'black'
                pos.color <- tri.colors[2]
            } else {
                neu.color <- tri.colors[2]
                pos.color <- tri.colors[3]
            }
            enrich.palette <- colorRampPalette(c(neg.color, neu.color, 
                                                 pos.color), 
                                               bias=color.distr, 
                                               interpolate='spline')
        } else {
            enrich.palette <- colorRampPalette(c('blue', 'black', 'yellow'), 
                                               bias=color.distr, 
                                               interpolate='spline')
        }
    } else {
        ## allow for different base color to be used in non-pair plots
        ## parse
        if(hm.color != "default") {
            base.color <- unlist(strsplit(hm.color, ':'))
            if(length(base.color)!=1 && length(base.color)!=2){
                ## raise error
                stop('You have specified more than two colors for your heatmap.\n\
                    Please provide only one or two colors.\n \
                    For example; "red" provides a white to red gradient, while \n\
                    "yellow:red" provides a yellow to red gradient')
            }
            enrich.palette <- colorRampPalette(c(base.color[1], base.color[2]))
        } else {
            enrich.palette <- colorRampPalette(c('snow', 'red2'))    
        }
    }

    hm_cols <- ncol(enrichList[[1]])

    # Adjust X-axis tick position. In a heatmap, X-axis is [0, 1].
    # Assume xticks$pos is from 0 to N(>0).
    xticks$pos <- xticks$pos / tail(xticks$pos, n=1)  # scale to the same size.

    # Define a function to calculate color breaks.
    ColorBreaks <- function(max.e, min.e, bam.pair, ncolor) {
    # Args:
    #   max.e: maximum enrichment value to be mapped to color.
    #   min.e: minimum enrichment value to be mapped to color.
    #   bam.pair: boolean tag for bam-pair.
    #   ncolor: number of colors to use.
    # Returns: vector of color breaks.

        # If bam-pair is used, create breaks for positives and negatives 
        # separately. If log2 ratios are all positive or negative, use only 
        # half of the color space.
        if(bam.pair) {
            max.e <- ifelse(max.e > 0, max.e, 1)
            min.e <- ifelse(min.e < 0, min.e, -1)
            c(seq(min.e, 0, length.out=ncolor / 2 + 1),
              seq(0, max.e, length.out=ncolor / 2 + 1)[-1])
        } else {
            seq(min.e, max.e, length.out=ncolor + 1)
        }
    }

    if(grepl(",", color.scale)) {
        scale.pair <- unlist(strsplit(color.scale, ","))
        scale.min <- as.numeric(scale.pair[1])
        scale.max <- as.numeric(scale.pair[2])
        if(scale.min >= scale.max) {
            warning("Color scale min value is >= max value.\n")
        }
        flood.pts <- c(scale.min, scale.max)
        brk.use <- ColorBreaks(scale.max, scale.min, bam.pair, ncolor)
    }

    # If color scale is global, calculate breaks and quantile here.
    if(color.scale == 'global') {
        flood.pts <- quantile(c(enrichList, recursive=T), c(flood.q, 1-flood.q))
        brk.use <- ColorBreaks(flood.pts[2], flood.pts[1], bam.pair, ncolor)
    }

    # Go through each unique region. 
    # Do NOT use "dopar" in the "foreach" loops here because this will disturb
    # the image order.
    go.list <- vector('list', length(uniq.reg))
    go.cluster <- vector('list', length(uniq.reg))

    names(go.list) <- uniq.reg
    names(go.cluster) <- uniq.reg

    for(ur in uniq.reg) {
        # ur <- uniq.reg[i]
        plist <- which(reg.list==ur)  # get indices in the config file.

        # Combine all profiles into one.
        # enrichCombined <- do.call('cbind', enrichList[plist])
        enrichSelected <- enrichList[plist]

        # If color scale is region, calculate breaks and quantile here.
        if(color.scale == 'region') {
            flood.pts <- quantile(c(enrichSelected, recursive=T), 
                                  c(flood.q, 1-flood.q))
            brk.use <- ColorBreaks(flood.pts[2], flood.pts[1], bam.pair, ncolor)
        }

        # Order genes.
        if(is.matrix(enrichSelected[[1]]) && nrow(enrichSelected[[1]]) > 1) {
            if(bam.pair) {
                lowCutoffs <- 0
            } else {
                lowCutoffs <- v.low.cutoff[plist]
            }
            g.order.list <- OrderGenesHeatmap(enrichSelected, lowCutoffs, 
                                              go.algo, go.paras)
            g.order <- g.order.list[[1]]
            g.cluster <- g.order.list[[2]]
            if(is.null(g.cluster)) {
                go.cluster[[ur]] <- NA
            } else{
                go.cluster[[ur]] <- rev(g.cluster[g.order])
            }
            go.list[[ur]] <- rev(rownames(enrichSelected[[1]][g.order, ]))
        } else {
            go.cluster[[ur]] <- NULL
            go.list[[ur]] <- NULL
        }

        if(!do.plot) {
            next
        }
  
        # Go through each sample and do plot.
        for(pj in plist) {
            if(!is.null(g.order)) {
                enrichList[[pj]] <- enrichList[[pj]][g.order, ]
            }

            # If color scale is local, calculate breaks and quantiles here.
            if(color.scale == 'local') {
                flood.pts <- quantile(c(enrichList[[pj]], recursive=T), 
                                      c(flood.q, 1-flood.q))
                brk.use <- ColorBreaks(flood.pts[2], flood.pts[1], bam.pair, 
                                       ncolor)
            }

            # Flooding extreme values.
            enrichList[[pj]][ enrichList[[pj]] < flood.pts[1] ] <- flood.pts[1]
            enrichList[[pj]][ enrichList[[pj]] > flood.pts[2] ] <- flood.pts[2]

            # Draw colorkey.
            image(z=matrix(brk.use, ncol=1), col=enrich.palette(ncolor), 
                  breaks=brk.use, axes=F, useRaster=T, main='Colorkey')
            axis(1, at=seq(0, 1, length.out=5), 
                 labels=format(brk.use[seq(1, ncolor + 1, length.out=5)], 
                               digits=1), 
                 lwd=1, lwd.ticks=1)

            # Draw heatmap.
            image(z=t(enrichList[[pj]]), col=enrich.palette(ncolor), 
                  breaks=brk.use, axes=F, useRaster=T, main=title2plot[pj])

            axis(1, at=xticks$pos, labels=xticks$lab, lwd=1, lwd.ticks=1)
        }
    }
    list(go.list,go.cluster)
}

trim <- function(x, p){
# Trim a numeric vector on both ends.
# Args:
#   x: numeric vector.
#   p: percentage of data to trim.
# Return: trimmed vector.

    low <- quantile(x, p)
    hig <- quantile(x, 1 - p)
    x[x > low & x < hig]
}

CalcSem <- function(x, rb=.05){ 
# Calculate standard error of mean for a numeric vector.
# Args:
#   x: numeric vector
#   rb: fraction of data to trim before calculating sem.
# Return: a scalar of the standard error of mean

    if(rb > 0){
        x <- trim(x, rb)
    }
    sem <- sd(x) / sqrt(length(x))
    ifelse(is.na(sem), 0, sem)
    # NOTE: this should be improved to handle exception that "sd" calculation 
    # emits errors.
}




## Leave for future reference.
#
# Set the antialiasing.
# type <- NULL
# if (capabilities()["aqua"]) {
#   type <- "quartz"
# } else if (capabilities()["cairo"]) {
#   type <- "cairo"
# } else if (capabilities()["X11"]) {
#   type <- "Xlib"
# }
# Set the output type based on capabilities.
# if (is.null(type)){
#   png(plot.name, width, height, pointsize=pointsize)

# } else {
#   png(plot.name, width, height, pointsize=pointsize, type=type)
# }
