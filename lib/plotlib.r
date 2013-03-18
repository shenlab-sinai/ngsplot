#### plotlib.r ####
# This contains the library for plotting related functions.
#
# Authors: Li Shen, Ningyi Shao
# 
# Created: Feb 19, 2013
# Last updated: Feb 21, 2013
#

SetupHeatmapDevice <- function(reg.list, uniq.reg, ng.list, pts, 
                               unit.width=4, reduce.ratio=30) {
# Configure parameters for heatmap output device. The output is used by 
# external procedures to setup pdf device ready for heatmap plotting.
# Args:
#   reg.list: region list as in config file.
#   uniq.reg: unique region list.
#   ng.list: number of genes per heatmap in the order as config file.
#   pts: data points (number of columns of heatmaps).
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

    # Setup image size.
    hm.width <- unit.width * max(reg.np)
    ipl <- .2 # inches per line. Obtained from par->'mai', 'mar'.
    m.bot <- 2; m.lef <- 1; m.top <- 2; m.rig <- 1 # margin size in lines.
    # Convert #gene to image height.
    reg.hei <- sapply(reg.ng, function(r) {
                    r * unit.width / pts / reduce.ratio + 
                    m.bot * ipl + m.top * ipl  # margins are like intercepts.
                })
    hm.height <- sum(reg.hei)

    # Setup layout of the heatmaps.
    lay.mat <- matrix(0, ncol=max(reg.np), nrow=length(reg.np))
    f.sta <- 1 # figure number start.
    for(i in 1:length(reg.np)){
        f.end <- f.sta + reg.np[i] - 1 # figure number end.
        lay.mat[i, 1:reg.np[i]] <- f.sta : f.end
        f.sta <- f.sta + reg.np[i]
    }

    list(reg.hei=reg.hei, hm.width=hm.width, hm.height=hm.height, 
         lay.mat=lay.mat, heatmap.mar=c(m.bot, m.lef, m.top, m.rig))
}

SetPtsSpline <- function(pint, lgint) {
# Set data points for spline function.
# Args:
#   pint: tag for point interval.
# Return: list of data points, middle data points, flanking data points.

    if(pint){  # point interval.
        pts <- 101  # data points to plot and to calculate SEM.
        m.pts <- 1  # middle part points.
        f.pts <- 50 # flanking part points.
    } else {
        pts <- 100
        if(lgint) {
            m.pts <- pts / 5 * 3
            f.pts <- pts / 5
        } else {
            m.pts <- pts / 5
            f.pts <- pts / 5 * 2
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

    regcovMat <- matrix(0, nrow=pts, ncol=nrow(ctg.tbl))
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
        confiMat <- matrix(0, nrow=pts, ncol=nrow(ctg.tbl))
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

genXticks <- function(reg2plot, pint, lgint, pts, flanksize, flankfactor){
# Generate X-ticks for plotting.
# Args:
#   reg2plot: string representation of region.
#   pint: point interval.
#   lgint: tag for large interval.
#   pts: data points.
#   flanksize: flanking region size in bps.
#   flankfactor: flanking region factor.
# Return: list of x-tick position and label.

    if(pint){   # point interval.
        lab.tbl <- c('TSS', 'TES', 'Center')
        names(lab.tbl) <- c('tss', 'tes', 'bed')
        mid.lab <- lab.tbl[reg2plot]
        tick.pos <- c(1, (pts - 1)/4, (pts - 1)/2 + 1, (pts - 1)/4*3, pts)
        tick.lab <- as.character(c(-flanksize, -flanksize/2, mid.lab, 
            flanksize/2, flanksize))
    }else{
        if(reg2plot == 'genebody'){
            left.lab <- 'TSS'
            right.lab <- 'TES'
        }else if(reg2plot == 'exon'){
            left.lab <- 'Acceptor'
            right.lab <- 'Donor'
        }else if(reg2plot == 'cgi' || reg2plot == 'bed'){
            left.lab <- 'Left'
            right.lab <- 'Right'
        }
        tick.pos <- c(0, pts/5, pts/5*2, pts/5*3, pts/5*4, pts)
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

plotmat <- function(regcovMat, title2plot, xticks, pts, m.pts, f.pts, pint,
                    shade.alp=0, confiMat=NULL, mw=1){
# Plot avg. profiles and standard errors around them.
# Args:
#   regcovMat: matrix for avg. profiles.
#   title2plot: profile names, will be shown in figure legend.
#   xticks: as is
#   pts: data points
#   m.pts: middle part data points
#   f.pts: flanking part data points
#   pint: tag for point interval
#   shade.alp: shading area alpha
#   confiMat: matrix for standard errors.

    # Smooth avg. profiles if specified.
    if(mw > 1){
        regcovMat <- as.matrix(runmean(regcovMat, k=mw, alg='C', 
                                       endrule='mean'))
    }

    # Choose colors.
    ncurve <- ncol(regcovMat)
    if(ncurve <= 8) {
        suppressMessages(require(RColorBrewer, warn.conflicts=F))
        col2use <- brewer.pal(ifelse(ncurve >= 3, ncurve, 3), 'Dark2')
        col2use <- col2use[1:ncurve]
    } else {
        col2use <- rainbow(ncurve)
    }
    col2use <- col2alpha(col2use, 0.8)

    # Plot profiles.
    ytext <- "Normalized Coverage (RPM)"
    xrange <- 1:pts
    matplot(xrange, regcovMat, xaxt='n', type="l", col=col2use, 
            lty="solid", lwd=3,
            xlab="Genomic Region (5' -> 3')", ylab=ytext)
    axis(1, at=xticks$pos, labels=xticks$lab, lwd=3, lwd.ticks=3)

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

    # If not point interval, add an extra vertical line.
    if(!pint){
        abline(v=f.pts, col="gray", lwd=2)
    }

    # Add gray lines indicating feature boundaries.
    abline(v=f.pts + m.pts, col="gray", lwd=2)

    # Legend.
    legend("topright", title2plot, text.col=col2use)
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

OrderGenesHeatmap <- function(n, enrichCombined, 
                    method=c('total', 'max', 'prod', 'diff', 'hc', 'pca')) {
# Order genes in combined heatmap data.
# Args: 
#   n: number of plots(such as histone marks) in the combined data.
#   enrichCombined: combined heatmap data.
#   method: algorithm used to order genes.
# Return: list of vectors of gene orders. In case of PCA, it may return more 
#   than one vectors of gene orders. Otherwise, the list length is 1.

    npts <- ncol(enrichCombined) / n  # number of data points for each profile.
    
    if(method == 'hc') {  # hierarchical clustering
        # Filter genes with zero sd.
        g.sd <- apply(enrichCombined, 1, sd)
        g.nz <- which(g.sd > 0)
        g.ze <- which(g.sd == 0)
        enrichCombined <- enrichCombined[g.nz, ]
        # Clustering and order genes.
        hc <- hclust(as.dist(1-cor(t(enrichCombined))), method='complete')
        # Notes: do NOT forget hc is applied to non-zero sd genes only.
        # The original gene indices must be recovered before return values.
        list(hc=c(g.nz[hc$order], g.ze))
    } else if(method == 'total') {  # overall enrichment of the 1st profile.
        list(total=order(rowSums(enrichCombined[, 1:npts])))
    } else if(method == 'max') {  # peak enrichment value of the 1st profile.
        list(max=order(apply(enrichCombined[, 1:npts], 1, max)))
    } else if(method == 'prod') {  # product of all profiles.
        g.prod <- foreach(r=iter(enrichCombined, by='row'), .combine='c', 
                            .multicombine=T, .maxcombine=1000) %dopar% {
            foreach(i=icount(n),  # go through each profile.
                    .combine='prod', .multicombine=T) %do% {
                col.sta <- (i - 1) * npts + 1
                col.end <- i * npts
                sum(r[col.sta:col.end], na.rm=T)
            }
        }
        list(prod=order(g.prod))
    } else if(method == 'diff') {  # difference between 1st and 2nd profiles.
        if(n < 2) {  # if only one profile available, return org. order.
            list(none=1:nrow(enrichCombined))
        } else {
            list(diff=order(rowSums(enrichCombined[, 1:npts]) - 
                rowSums(enrichCombined[, (npts + 1):(npts * 2)])))
        }
    } else if(method == 'pca') {  # principal component analysis.
        # Reduce the data to a small number of bins per profile.
        nbin <- 10
        enrich.reduced <- foreach(i=icount(n), .combine='cbind', 
                                .multicombine=T) %do% {
            # Go through each profile.
            col.sta <- (i - 1) * npts + 1
            col.end <- i * npts
            # Column breaks represent bin boundaries.
            col.breaks <- seq(col.sta, col.end, length.out=nbin + 1)
            foreach(j=icount(nbin), .combine='cbind', .multicombine=T) %dopar% {
                # Go through each bin.
                rowSums(enrichCombined[, col.breaks[j]:col.breaks[j + 1]])
            }
        }
        # Pull out all pc's that equal at least 10% variance of the 1st pc.
        enrich.pca <- prcomp(enrich.reduced, center=F, scale=F, tol=sqrt(.1))
        # Order genes according to each pc.
        pc.order <- foreach(i=icount(ncol(enrich.pca$x))) %dopar% {
            order(enrich.pca$x[, i])
        }
        names(pc.order) <- paste('pc', 1:ncol(enrich.pca$x), sep='')
        pc.order
    }
}


plotheat <- function(reg.list, uniq.reg, enrichList, go.algo, title2plot, 
                     xticks, flood.q=.02) {
# Plot heatmaps with genes ordered according to some algorithm.
# Args:
#   reg.list: factor vector of regions as in configuration.
#   uniq.reg: character vector of unique regions.
#   enrichList: list of heatmap data.

    # Setup basic parameters.
    ncolor <- 256
    enrich.palette <- colorRampPalette(c('snow', 'red2'))
    hm_cols <- ncol(enrichList[[1]])

    # Adjust X-axis tick position. In a heatmap, X-axis is [0, 1].
    mat.xsz <- tail(xticks$pos, n=1) - xticks$pos[1] + 1
    xticks$pos <- xticks$pos - xticks$pos[1] + 1  # shift left-most to 1.
    xticks$pos <- xticks$pos / mat.xsz  # scale to the same size.

    # Go through each unique region. 
    # Do NOT use "dopar" in the "foreach" loops here because this will disturb
    # the image order.
    for(i in 1:length(uniq.reg)) {
        ur <- uniq.reg[i]
        plist <- which(reg.list==ur)    # get indices in the config file.

        # Combine all profiles into one.
        enrichCombined <- do.call('cbind', enrichList[plist])

        # Order genes.
        if(go.algo != 'none' && nrow(enrichCombined) > 1) {
            g.order <- OrderGenesHeatmap(length(plist), enrichCombined, go.algo)
            enrichCombined <- enrichCombined[g.order[[1]], ]
        }
        # for now, just use the 1st gene order.
  
        # Split combined profiles back into individual heatmaps.
        for(j in 1:length(plist)) {
            pj <- plist[j]  # index in the original config.

            enrichList[[pj]] <- enrichCombined[, ((j-1)*hm_cols+1) : 
                                                (j*hm_cols)]

            # Flooding extreme values which are identified by quantiles.
            flood.pts <- quantile(c(enrichList[[pj]], recursive=T), 
                                c(flood.q, 1-flood.q))
            enrichList[[pj]][ enrichList[[pj]] < flood.pts[1] ] <- flood.pts[1]
            enrichList[[pj]][ enrichList[[pj]] > flood.pts[2] ] <- flood.pts[2]

            # Draw heatmap.
            image(z=t(enrichList[[pj]]), col=enrich.palette(ncolor), axes=F, 
                  useRaster=T, main=title2plot[pj])

            axis(1, at=xticks$pos, labels=xticks$lab, lwd=1, lwd.ticks=1)
        }
    }
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
