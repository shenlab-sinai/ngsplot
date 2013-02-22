# Function to check if the range exceeds coverage vector boundaries.
checkBound <- function(start, end, range, chrlen){
    if(end + range > chrlen || start - range < 1)
        return(FALSE)   # out of boundary.
    else
        return(TRUE)
}

Spline3Sec <- function(cov.list, pts.list) {
# For 3 sections of continuous coverage, spline each according to specified 
# data points and return a concatenated, interpolated curve.
# Args:
#   cov.list: a list of three coverage sections.
#   pts.list: a list of three integers for data points.
# Notes: the names of cov.list and pts.list must be the same for index purpose.
# Suggested: use 'l', 'm', and 'r' for the names.

    start <- 1
    end <- length(cov.list[[1]])

    xp <- foreach(s=names(cov.list[-1]),  # start from the 2nd cov in the list.
            .init=seq(start, end, length.out=pts.list[[1]]), 
            .combine='c') %do% {
        start <- end + 1
        end <- end + length(cov.list[[s]])
        # Must use "ceiling" to prevent return 0.
        ceiling(seq(start, end, length.out=pts.list[[s]]))
    }

    all.cov <- do.call('c', cov.list)

    spline(1:length(all.cov), all.cov, xout=xp)$y
}

# Extract and interpolate coverage vector from a genomic region with 3 sections:
# Upstream flanking, middle and downstream flanking.
extrCov3Sec <- function(chrcov, start, end, flanking, m.pts, f.pts, strand){

    left.cov <- as.vector(seqselect(chrcov, start - flanking, start - 1))
    right.cov <- as.vector(seqselect(chrcov, end + 1, end + flanking))
    middle.cov <- as.vector(seqselect(chrcov, start, end))

    # # Sample middle and flanking regions separately.
    # left.cov <- spline(1:length(left.cov), left.cov, n=f.pts)$y
    # right.cov <- spline(1:length(right.cov), right.cov, n=f.pts)$y
    # middle.cov <- spline(1:length(middle.cov), middle.cov, n=m.pts)$y

    # all.cov <- c(left.cov, middle.cov, right.cov)

    splined.cov <- Spline3Sec(list(l=left.cov, m=middle.cov, r=right.cov),
                              list(l=f.pts, m=m.pts, r=f.pts))

    if(strand == '+') {
        splined.cov
    } else {
        rev(splined.cov)
    }
}

# Extract and interpolate coverage vector from a genomic region.
extrCovSec <- function(chrcov, start, end, flanking, pts, strand, weight){
    cov <- as.vector(seqselect(chrcov, start - flanking, end + flanking))
    cov.intp <- spline(1:length(cov), cov, n=pts)$y
    if(weight){
        cov.intp <- (length(cov) / pts) * cov.intp
    }
    if(strand == '+'){
        cov.intp
    }else{
        rev(cov.intp)
    }
}


# Extract and concatenate coverages for a mRNA using exon model. Then do interpolation.
extrCovExons <- function(chrcov, start, end, ranges, flanking, m.pts, f.pts, 
                        strand){

    exon.cov <- as.vector(seqselect(chrcov, ranges))
    left.cov <- as.vector(seqselect(chrcov, start - flanking, start - 1))
    right.cov <- as.vector(seqselect(chrcov, end + 1, end + flanking))

    splined.cov <- Spline3Sec(list(l=left.cov, m=exon.cov, r=right.cov),
                              list(l=f.pts, m=m.pts, r=f.pts))

    if(strand == '+') {
        splined.cov
    } else {
        rev(splined.cov)
    }

    # # Sample from spline function.
    # exon.cov <- spline(1:length(exon.cov), exon.cov, n=m.pts)$y
    # left.cov <- spline(1:length(left.cov), left.cov, n=f.pts)$y
    # right.cov <- spline(1:length(right.cov), right.cov, n=f.pts)$y

    # all.cov <- c(left.cov, exon.cov, right.cov)

    # if(strand == '+'){
    #     all.cov
    # }else{
    #     rev(all.cov)
    # }
}

# Extract coverage vector from a genomic region with a middle point and symmetric flanking regions.
extrCovMidp <- function(chrcov, midp, flanking, pts, strand){
    res.cov <- as.vector(seqselect(chrcov, midp - flanking, midp + flanking))
    res.cov <- spline(1:length(res.cov), res.cov, n=pts)$y
    if(strand == '+'){
        res.cov
    }else{
        rev(res.cov)
    }
}

do.par.cov <- function(rec, read.coverage.n, flanksize, flankfactor, 
                        m.pts, f.pts, reg2plot, pint, exon.ranges=NULL){

    chrom <- as.character(rec$chrom)
    if(!chrom %in% names(read.coverage.n)){
        return(NULL)
    }
    strand <- rec$strand

    if(!is.null(exon.ranges)){
        start <- start(exon.ranges[1])
        end <- end(tail(exon.ranges, n=1))
        reglen <- sum(end(exon.ranges) - start(exon.ranges) + 1)
    } else {
        start <- rec$start
        end <- rec$end
        reglen <- end - start + 1
    }

    if(flankfactor > 0 && !pint){
        flanksize <- round(reglen * flankfactor)
    }

    chr.cov <- read.coverage.n[[chrom]]

    if((reg2plot == 'tss' && strand == '+') || 
        (reg2plot == 'tes' && strand == '-')){
        
        if(!checkBound(start, start, flanksize, length(chr.cov)))
            return(NULL)

        result <- extrCovMidp(chr.cov, start, flanksize, m.pts + f.pts*2, strand)

    } else if(reg2plot == 'tss' && strand == '-' || 
            reg2plot == 'tes' && strand == '+'){

        if(!checkBound(end, end, flanksize, length(chr.cov)))
            return(NULL)

        result <- extrCovMidp(chr.cov, end, flanksize, m.pts + f.pts*2, strand)

    } else{
        if(!checkBound(start, end, flanksize, length(chr.cov)))
            return(NULL)

        if(!is.null(exon.ranges)){  # RNA-seq plot using exon model.
            result <- extrCovExons(chr.cov, start, end, exon.ranges, flanksize, 
                                    m.pts, f.pts, strand)
        } else{
            result <- extrCov3Sec(chr.cov, start, end, flanksize, m.pts, f.pts, 
                                strand)
        }
    }

    result
}
