# Create semi-transparent colors.
col2alpha <- function(col2use, alpha){
	apply(col2rgb(col2use), 2, function(x){
		rgb(x[1], x[2], x[3], 
			alpha=alpha*255, 
			maxColorValue=255)
		})
}

# Smooth plotting data in a vector.
smoothvec <- function(v, radius, method=c('mean', 'median')){
	stopifnot(is.vector(v))
	stopifnot(length(v) > 0)
	stopifnot(radius > 0 && radius < 1)
	halfwin <- ceiling(length(v) * radius)
	s <- vector()
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

# Smooth plotting data which can be a matrix or a vector.
smoothplot <- function(m, radius, method=c('mean', 'median')){
	stopifnot(is.matrix(m) || is.vector(m))
	if(is.matrix(m)){
		for(i in 1:ncol(m)){
			m[, i] <- smoothvec(m[, i], radius, method)
		}
	}else{
		m <- smoothvec(m, radius, method)
	}
	m
}


# Generate X-axis ticks.
gen_xticks <- function(reg2plot, intsize, flanksize, flankfactor){
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
	list(pos=tick.pos, lab=tick.lab)
}

# Function to plot average profiles.
plotmat <- function(regcovMat, title2plot, xticks, shade.alp=0, confiMat=NULL){
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
	axis(1, at=xticks$pos, labels=xticks$lab, lwd=3, lwd.ticks=3)

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

	# Render standard errors around each curve.
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

	# Add boundaries and legend.
	abline(v=-(intsize-1)/2, col="gray", lwd=2)
	abline(v=(intsize-1)/2, col="gray", lwd=2)
	legend("topright", title2plot, text.col=col2use)
}

# Use spline to sample each row of a matrix. Use this function for heatmap.
spline_mat <- function(mat, n=100){

	foreach(r=iter(mat, by='row'), 
			.combine='rbind', .multicombine=T) %dopar% {
		spline(1:length(r), r, n)$y
	}
}

# Function to plot heatmaps.
plotheat <- function(reg.list, uniq.reg, enrichList, title2plot, xticks,
						flood.q=.01){
	# Setup basic parameters.
    ncolor <- 256
    enrich.palette <- colorRampPalette(c('snow', 'red2'))
    hm_cols <- ncol(enrichList[[1]])

    # Adjust X-axis tick position. In a heatmap, X-axis is [0, 1].
    mat.xsz <- tail(xticks$pos, n=1) - xticks$pos[1] + 1
    xticks$pos <- xticks$pos - xticks$pos[1] + 1 # shift left-most to 1.
    xticks$pos <- xticks$pos / mat.xsz # scale to the same size.

    # Go through each unique region. 
    # Do NOT use "dopar" in the "foreach" loops here because this will disturb
    # the image order.
    foreach(ur=iter(uniq.reg)) %do% {

    	plist <- which(reg.list==ur)	# get indices in the config file.
	    # Combine all profiles into one.
	    enrichCombined <- do.call('cbind', enrichList[plist])

	    # Filter genes with zero sd.
	    g.sd <- apply(enrichCombined, 1, sd)
	    enrichCombined <- enrichCombined[g.sd > 0, ]
	    
	    # Clustering and order genes.
	    hc <- hclust(as.dist(1-cor(t(enrichCombined))), method='complete')
	    enrichCombined <- enrichCombined[hc$order, ]
	    
	    # Split combined profiles back into individual heatmaps.
    	foreach(j=1:length(plist)) %do% {
    		pj <- plist[j]	# index in the original config.

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


## Leave for future reference.
#
# Set the antialiasing.
# type <- NULL
# if (capabilities()["aqua"]) {
# 	type <- "quartz"
# } else if (capabilities()["cairo"]) {
# 	type <- "cairo"
# } else if (capabilities()["X11"]) {
# 	type <- "Xlib"
# }
# Set the output type based on capabilities.
# if (is.null(type)){
# 	png(plot.name, width, height, pointsize=pointsize)

# } else {
# 	png(plot.name, width, height, pointsize=pointsize, type=type)
# }
