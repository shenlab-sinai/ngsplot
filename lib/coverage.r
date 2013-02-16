# Function to check if the range exceeds coverage vector boundaries.
checkBound <- function(start, end, range, chrlen){
    if(end + range > chrlen || start - range < 1)
        return(FALSE)   # out of boundary.
    else
        return(TRUE)
}

# Extract and interpolate coverage vector from a genomic region with 3 sections:
# 5' raw region, variable middle region and 3' raw region.
extrCov3Sec <- function(chrcov, start, end, ninterp, flanking, strand, weight){
    left.cov <- as.vector(seqselect(chrcov, start - flanking, start - 1))
    right.cov <- as.vector(seqselect(chrcov, end + 1, end + flanking))
    middle.cov <- as.vector(seqselect(chrcov, start, end))
    middle.cov.intp <- spline(1:length(middle.cov), middle.cov, n=ninterp)$y
    if(weight){
        middle.cov.intp <- (length(middle.cov) / ninterp) * middle.cov.intp
    }
    if(strand == '+'){
        return(c(left.cov, middle.cov.intp, right.cov))
    }else{
        return(rev(c(left.cov, middle.cov.intp, right.cov)))
    }
}

# Extract and interpolate coverage vector from a genomic region.
extrCovSec <- function(chrcov, start, end, ninterp, flanking, strand, weight){
    cov <- as.vector(seqselect(chrcov, start - flanking, end + flanking))
    cov.intp <- spline(1:length(cov), cov, n=ninterp)$y
    if(weight){
        cov.intp <- (length(cov) / ninterp) * cov.intp
    }
    if(strand == '+'){
        return(cov.intp)
    }else{
        return(rev(cov.intp))
    }
}


# Extract and concatenate coverages for a mRNA using exon model. Then do interpolation.
extrCovExons <- function(chrcov, ranges, ninterp, strand, weight){
    cov <- as.vector(seqselect(chrcov, ranges))
    cov.intp <- spline(1:length(cov), cov, n=ninterp)$y
    if(weight){
        cov.intp <- (length(cov) / ninterp) * cov.intp
    }
    if(strand == '+'){
        return(cov.intp)
    }else{
        return(rev(cov.intp))
    }
}

# Extract coverage vector from a genomic region with a middle point and symmetric flanking regions.
extrCovMidp <- function(chrcov, midp, flanking, strand){
    res.cov <- as.vector(seqselect(chrcov, midp - flanking, midp + flanking))
    if(strand == '+'){
        return(res.cov)
    }else{
        return(rev(res.cov))
    }
}

do.par.cov <- function(k, plot.coord, read.coverage.n, rnaseq.gb,
                     flankfactor, reg2plot, genemodel, weight.genlen,
                     intsize, old_flanksize, flanksize){
    chrom <- as.character(plot.coord[k, ]$chrom)
    if(!chrom %in% names(read.coverage.n)) return(NULL)
    strand <- plot.coord[k, ]$strand
    if(flankfactor > 0 && !rnaseq.gb){
        flanksize <- floor((plot.coord[k, ]$end - plot.coord[k, ]$start + 1)*flankfactor)
    }
    if((reg2plot == 'tss' && strand == '+') || (reg2plot == 'tes' && strand == '-')){
        if(!checkBound(plot.coord[k, ]$start, plot.coord[k, ]$start, flanksize, length(read.coverage.n[[chrom]])))
            return(NULL)
        result <- extrCovMidp(read.coverage.n[[chrom]], plot.coord[k, ]$start, flanksize, strand)
    }else if(reg2plot == 'tss' && strand == '-' || reg2plot == 'tes' && strand == '+'){
        if(!checkBound(plot.coord[k, ]$end, plot.coord[k, ]$end, flanksize, length(read.coverage.n[[chrom]])))
            return(NULL)
        result <- extrCovMidp(read.coverage.n[[chrom]], plot.coord[k, ]$end, flanksize, strand)
    }else{
        if(!checkBound(plot.coord[k, ]$start, plot.coord[k, ]$end, flanksize, length(read.coverage.n[[chrom]])))
            return(NULL)
        if(rnaseq.gb){  # RNA-seq plot using exon model.
            exon.ranges <- genemodel$exonmodel[[plot.coord$tid[k]]]$ranges
            result <- extrCovExons(read.coverage.n[[chrom]], exon.ranges, intsize, strand, weight.genlen)
        }else{
            if(flankfactor > 0){    # one section coverage.
                result <- extrCovSec(read.coverage.n[[chrom]], plot.coord[k, ]$start, plot.coord[k, ]$end, intsize+2*old_flanksize, flanksize, strand, weight.genlen)
            }else{  # three section coverage.
                result <- extrCov3Sec(read.coverage.n[[chrom]], plot.coord[k, ]$start, plot.coord[k, ]$end, intsize, flanksize, strand, weight.genlen)
            }
        }
    }
    result
}
