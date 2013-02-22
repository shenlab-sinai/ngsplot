# Create semi-transparent colors.
col2alpha <- function(col2use, alpha){
	 col2use <- apply(col2rgb(col2use), 2, function(x){rgb(x[1], x[2], x[3], alpha=alpha*255, maxColorValue=255)})
	 col2use
}


# Function to plot average profiles.
plotmat <- function(plot.name, width, height, pointsize, 
	reg2plot, flanksize, intsize, flankfactor, shade.alp, rnaseq.gb,
	regcovMat, title2plot, confiMat=NULL){
	ppi <- 256

	# Plot into image file.
	# Set the antialiasing.
	type <- NULL
	if (capabilities()["aqua"]) {
		type <- "quartz"
	} else if (capabilities()["cairo"]) {
		type <- "cairo"
	} else if (capabilities()["X11"]) {
		type <- "Xlib"
	}
	pdf(plot.name, width=width/ppi, height=height/ppi)

	# Set the output type based on capabilities.
	# if (is.null(type)){
	# 	png(plot.name, width, height, pointsize=pointsize)

	# } else {
	# 	png(plot.name, width, height, pointsize=pointsize, type=type)
	# }

	# Choose colors.
	ncurve <- ncol(regcovMat)
	if(ncurve <= 8) {
		library(RColorBrewer)
		col2use <- brewer.pal(ifelse(ncurve >= 3, ncurve, 3), 'Dark2')
		col2use <- col2use[1:ncurve]
	} else {
		col2use <- rainbow(ncurve)
	}
	col2use <- col2alpha(col2use, 0.8)
	# Draw curves.
	ytext <- "Normalized Coverage(RPM)"
	xrange <- ((-flanksize-(intsize-1)/2) : (flanksize+(intsize-1)/2))
	matplot(xrange, regcovMat, xaxt='n', type="l", col=col2use, lty="solid", lwd=3,
		xlab='DNA basepair(or interpolated)', ylab=ytext)
	# Handle ticks.
	if(reg2plot == 'tss' || reg2plot == 'tes'){
		mid.lab <- ifelse(reg2plot == 'tss', 'TSS', 'TES')
		tick.pos <- c(-flanksize, -flanksize/2, 0, flanksize/2, flanksize)
		tick.lab <- as.character(c(-flanksize, -flanksize/2, mid.lab, flanksize/2, flanksize))
	}else{
		if(reg2plot == 'genebody'){
			left.lab <- 'TSS'; right.lab <- 'TES'
		}else if(reg2plot == 'exon'){
			left.lab <- 'Acceptor'; right.lab <- 'Donor'
		}else if(reg2plot == 'cgi' || reg2plot == 'bed'){
			left.lab <- 'Left'; right.lab <- 'Right'
		}
		if(rnaseq.gb){
			tick.pos <- c(-(intsize-1)/2, (intsize-1)/2)
			tick.lab <- c(left.lab, right.lab)
		}else{
			tick.pos <- c((-flanksize-(intsize-1)/2), (-flanksize/2-(intsize-1)/2),
				-(intsize-1)/2, (intsize-1)/2,
				(flanksize/2+(intsize-1)/2), (flanksize+(intsize-1)/2))
			if(flankfactor > 0){
				tick.lab <- c(paste(c(-floor(flankfactor*100), -floor(flankfactor*50)), '%', sep=''),
				left.lab, right.lab,
				paste(c(floor(flankfactor*50), floor(flankfactor*100)), '%', sep=''))
			}else{
				tick.lab <- as.character(c(-flanksize, -flanksize/2,
					left.lab, right.lab,
					flanksize/2, flanksize))
			}
		}
	}
	axis(1, at=tick.pos, labels=tick.lab, lwd=3, lwd.ticks=3)
	# Add shaded area.
	if(shade.alp > 0){
		for(i in 1:ncol(regcovMat)){
			v.x <- c(xrange[1], xrange, xrange[length(xrange)])
			v.y <- regcovMat[, i]
			v.y <- c(0, v.y, 0)
			col.rgb <- col2rgb(col2use[i])
			p.col <- rgb(col.rgb[1,1], col.rgb[2,1], col.rgb[3,1], alpha=shade.alp*255, maxColorValue=255)
			polygon(v.x, v.y, density=-1, border=NA, col=p.col)
		}
	}
	if(!is.null(confiMat)){
		for(i in 1:ncol(confiMat)){
			step <- length(regcovMat[,i])/length(confiMat[,i])
			pos <- xrange[round((0:length(confiMat[,i]))*step)]
			pos.y <- round((0:length(confiMat[,i]))*step)
			v.x <- c(pos, rev(pos))
			v.y <- c((regcovMat[pos.y,i] + confiMat[,i]), rev(regcovMat[pos.y,i] - confiMat[,i]))
			col.rgb <- col2rgb(col2use[i])
			p.col <- rgb(col.rgb[1,1], col.rgb[2,1], col.rgb[3,1], alpha=0.2*255, maxColorValue=255)
			polygon(v.x, v.y, density=-1, border=NA, col=p.col)
		}
	}
	abline(v=-(intsize-1)/2, col="gray", lwd=2)
	abline(v=(intsize-1)/2, col="gray", lwd=2)
	legend("topright", title2plot, text.col=col2use)
	dev.off()
}


# Function to plot heatmaps.
plotheat <- function(plot.name, width, height, pointsize, 
	enrichList, title2plot, flood.q=.01){

	# Setup basic parameters.
    ncolor <- 256
    enrich.palette <- colorRampPalette(c('white', 'red'))
    np <- length(enrichList)

    # Filter genes with zero sd.
    for(i in 1:np){
	    g.sd <- apply(enrichList[[i]], 1, sd)
	    enrichList[[i]] <- enrichList[[i]][g.sd > 0, ]
    }
    
    # Clustering and order genes.
    for(i in 1:np){
	    hc <- hclust(as.dist(1-cor(t(enrichList[[i]]))), method='complete')
	    enrichList[[i]] <- enrichList[[i]][hc$order, ]
	}
    
    # Flooding extreme values which are identified by quantiles.
    flood.pts <- quantile(c(mat, recursive=T), c(flood.q, 1-flood.q))
    mat[mat < flood.pts[1]] <- flood.pts[1]
    mat[mat > flood.pts[2]] <- flood.pts[2]

    # Create multi-panel image area and plot.
    par(mar=c(0.1, 0.1, 2, 0.1), mfcol=c(1, np))
    for(i in 1:np){
	    image(z=t(enrichList[[i]]), col=enrich.palette(ncolor), axes=F, 
	    	useRaster=T, main=title2plot[i])
	}

}