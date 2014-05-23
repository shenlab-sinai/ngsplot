# Function to check if the range exceeds coverage vector boundaries.
checkBound <- function(start, end, range, chrlen){
    if(end + range > chrlen || start - range < 1)
        return(FALSE)   # out of boundary.
    else
        return(TRUE)
}

Bin <- function(v, n) {
# Function to bin a coverage vector.
# Args:
#   v: coverage vector.
#   n: number of bins.
# Return: vector of binned values.

    bin.bound <- seq(1, length(v), length.out=n + 1)
    sapply(1:n, function(i) {
            mean(v[bin.bound[i]:bin.bound[i + 1]])
        })
}

SplineRev3Sec <- function(cov.list, v.fls, pts.list, v.strand, algo="spline") {
# For 3 sections of continuous coverage, spline each according to specified 
# data points and return concatenated, interpolated curves.
# Args:
#   cov.list: a list of coverage vectors. Each vector represents a gene.
#   v.fls: vector of flanking region size.
#   pts.list: a list of three integers for data points.
#   v.strand: factor vector of strands.
#   algo: algorithm used to normalize coverage vectors (spline, bin).
# Return: matrix of interpolated coverage: each row represents a gene; each 
#         column represents a data point. Coverage from exons are concatenated.
# Notes: the names of cov.list and pts.list must be the same for index purpose.
# Suggested: use 'l', 'm', and 'r' for the names.

    # Create an empty coveage matrix first.
    tot.pts <- sum(unlist(pts.list)) - 2
    cov.mat <- matrix(nrow=length(cov.list), ncol=tot.pts)

    for(i in 1:length(cov.list)) {
        left.cov <- head(cov.list[[i]], n=v.fls[i])
        right.cov <- tail(cov.list[[i]], n=v.fls[i])
        mid.cov <- window(cov.list[[i]], start=v.fls[i] + 1, 
                          width=length(cov.list[[i]]) - 2*v.fls[i])
        if(algo == "spline") {
            left.cov <- spline(1:length(left.cov), left.cov, n=pts.list$l)$y
            right.cov <- spline(1:length(right.cov), right.cov, n=pts.list$r)$y
            mid.cov <- spline(1:length(mid.cov), mid.cov, n=pts.list$m)$y
        } else if(algo == "bin") {
            left.cov <- Bin(left.cov, n=pts.list$l)
            right.cov <- Bin(right.cov, n=pts.list$r)
            mid.cov <- Bin(mid.cov, n=pts.list$m)
        } else {
            # pass.
        }
        # Fuse the two points at the boundary and concatenate.
        f1 <- (tail(left.cov, n=1) + head(mid.cov, n=1)) / 2
        f2 <- (tail(mid.cov, n=1) + head(right.cov, n=1)) / 2
        con.cov <- c(left.cov[-length(left.cov)], f1, 
                     mid.cov[-c(1, length(mid.cov))], 
                     f2, right.cov[-1])
        # browser()
        if(v.strand[i] == '+') {
            cov.mat[i, ] <- con.cov
        } else {
            cov.mat[i, ] <- rev(con.cov)
        }
    }

    cov.mat
}

extrCov3Sec <- function(v.chrom, v.start, v.end, v.fls, v.strand, m.pts, f.pts, 
                        bufsize, algo, ...) {
# Extract and interpolate coverage vectors from genomic regions with 3 sections.
# Args:
#   v.chrom: factor vector of chromosome names.
#   v.start: integer vector of region start.
#   v.end: integer vector of region end.
#   v.fls: integer vector of flanking region size in bps.
#   v.strand: factor vector of gene strands.
#   m.pts: data points for middle interval.
#   f.pts: data points for flanking region.
#   bufsize: integer; buffers are added to both ends of each region.
#   algo: algorithm used to normalize coverage vectors.
# Return: matrix of interpolated coverage: each row represents a gene; each 
#         column represents a data point.

    interflank.gr <- GRanges(seqnames=v.chrom, ranges=IRanges(
                             start=v.start - v.fls - bufsize, 
                             end=v.end + v.fls + bufsize))
    # interflank.cov <- covBam(bam.file, sn.inbam, interflank.gr, fraglen, 
    #                          map.qual, bowtie)
    interflank.cov <- covBam(interflank.gr, v.strand, ...)

    # Trim buffers from coverage vectors.
    interflank.cov <- lapply(interflank.cov, function(v) {
            window(v, start=bufsize + 1, width=length(v) - 2 * bufsize)
        })

    # Interpolate and reverse coverage vectors.
    SplineRev3Sec(interflank.cov, v.fls, list(l=f.pts, m=m.pts, r=f.pts), 
                  v.strand, algo)
}

sub3CovList <- function(all.cov, v.left, v.right) {
# For a list of coverage vectors, separate them into three lists.
# Args:
#   all.cov: list of coverage vectors.
#   v.left: integer vector of left flanking size.
#   v.right: integer vector of right flanking size.
# Return: list of lists of coverage vectors. The outer list is named as:
#         "l", "m", "r".

    left.cov <- foreach(i=icount(length(all.cov))) %do% {
        head(all.cov[[i]], n=v.left[i])
    }
    mid.cov <- foreach(i=icount(length(all.cov))) %do% {
        len <- length(all.cov[[i]])
        all.cov[[i]][ (v.left[i] + 1) : (len - v.right[i]) ]
    }
    right.cov <- foreach(i=icount(length(all.cov))) %do% {
        tail(all.cov[[i]], n=v.right[i])
    }
    list(l=left.cov, m=mid.cov, r=right.cov)
}

extrCovExons <- function(v.chrom, exonranges.list, v.fls, v.strand, 
                         m.pts, f.pts, bufsize, algo, ...) {
# Extract coverage vectors for transcripts with exon models.
# Args:
#   v.chrom: factor vector of chromosome names.
#   exonranges.list: list of IRanges objects for exon coordinates.
#   v.fls: integer vector of flanking region size in bps.
#   v.strand: factor vector of gene strands.
#   m.pts: data points for middle interval.
#   f.pts: data points for flanking region.
#   bufsize: integer; buffers are added to both ends of each region.
#   algo: algorithm used to normalize coverage vectors.
# Return: matrix of interpolated coverage: each row represents a gene; each 
#         column represents a data point. Coverage from exons are concatenated.

    # Construct ranges including exon and flanking regions.
    exonflank.list <- vector('list', length=length(exonranges.list))
    for(i in 1:length(exonranges.list)) {
        r.mod <- exonranges.list[[i]]  # IRanges object to be modified.
        n <- length(r.mod)
        # Add flanking and buffer regions.
        start(r.mod)[1] <- start(r.mod)[1] - v.fls[i] - bufsize
        end(r.mod)[n] <- end(r.mod)[n] + v.fls[i] + bufsize
        exonflank.list[[i]] <- GRanges(seqnames=v.chrom[i], ranges=r.mod)
    }

    exonflank.cov <- covBamExons(exonflank.list, v.strand, ...)

    # Trim buffers from coverage vectors.
    exonflank.cov <- lapply(exonflank.cov, function(v) {
            window(v, start=bufsize + 1, width=length(v) - 2 * bufsize)
        })

    SplineRev3Sec(exonflank.cov, v.fls, list(l=f.pts, m=m.pts, r=f.pts), 
                  v.strand, algo)
}


extrCovMidp <- function(v.chrom, v.midp, flanksize, v.strand, pts, bufsize, 
                        algo, ...) {
# Extract coverage vectors with a middle point and symmetric flanking regions.
# Args:
#   v.chrom: factor vector of chromosome names.
#   v.midp: integer vector of middle points.
#   flanksize: integer of flanking region size in bps.
#   v.strand: factor vector of gene strands.
#   pts: data points to spline into.
#   bufsize: integer; buffers are added to both ends of each region.
#   algo: algorithm used to normalize coverage vectors.
# Return: matrix of interpolated coverage: each row represents a gene; each 
#         column represents a data point.
    
    granges <- GRanges(seqnames=v.chrom, ranges=IRanges(
                       start=v.midp - flanksize - bufsize, 
                       end=v.midp + flanksize + bufsize))
    # cov.list <- covBam(bam.file, sn.inbam, granges, fraglen, map.qual, bowtie)
    # browser()
    cov.list <- covBam(granges, v.strand, ...)

    # Trim buffers from coverage vectors.
    cov.list <- lapply(cov.list, function(v) {
            window(v, start=bufsize + 1, width=length(v) - 2 * bufsize)
        })

    # Interpolate and reverse coverage vectors and assemble into a matrix.
    cov.mat <- matrix(nrow=length(cov.list), ncol=pts)
    # cov.mat <- matrix(nrow=length(cov.list), ncol=length(cov.list[[1]]))
    for(i in 1:length(cov.list)) {
        if(is.null(cov.list[[i]])) {
            cov.mat[i, ] <- vector('integer', length=pts)
            # cov.mat[i, ] <- vector('integer', length=length(cov.list[[1]]))
        } else {
            if(algo == "spline") {
                cov.mat[i, ] <- spline(1:length(cov.list[[i]]), cov.list[[i]], 
                                       n=pts)$y
            } else if(algo == "bin") {
                cov.mat[i, ] <- Bin(cov.list[[i]], n=pts)
            } else {
                # pass.
            }
            # cov.mat[i, ] <- cov.list[[i]]
            if(v.strand[i] == '-') {
                cov.mat[i, ] <- rev(cov.mat[i, ])
            }
        }
    }
    # browser()
    cov.mat
}

scanBamRevOrder <- function(org.gr, sbp) {
# ScanBamParam re-arranges the input genomic ranges. Use range info to 
# construct a string vector to find the order to reverse it.
    org.grnames <- paste(org.gr$seqnames, org.gr$start, org.gr$end, sep=':')
    sbw.gr <- as.data.frame(bamWhich(sbp))  # scan-bam-ed
    sbw.grnames <- paste(sbw.gr$space, sbw.gr$start, sbw.gr$end, sep=':')
    match(org.grnames, sbw.grnames)
}

genZeroList <- function(llen, v.vlen) {
# Generate a list of specific length with each element being an Rle vector of 
# zeros with specific length.
# Args:
#   llen: list length
#   v.vlen: vector of vector lengths.

    llen <- as.integer(llen)
    stopifnot(llen > 0)
    res <- vector('list', length=llen)
    for(i in 1:llen) {
        res[[i]] <- Rle(0, v.vlen[i])
    }
    res
}

covBam <- function(granges, v.strand, bam.file, sn.inbam, fraglen, 
                   map.qual=20, bowtie=F, extr.only=F, 
                   strand.spec=c('both', 'same', 'opposite')) {
# Extract coverage vectors from bam file for a list of genes.
# Args:
#   granges: list of GRanges objects representing genomic coordinates 
#            to extract coverge from. In the case of RNA-seq, each GRanges
#            object has multiple ranges representing exons and flanking 
#            regions. In the case of ChIP-seq, each GRanges object has one
#            range representing left boundary to right boundary.
#   v.strand: vector of strand info.
#   bam.file: character string refers to the path of a bam file.
#   sn.inbam: vector of chromosome names in the bam file.
#   fraglen: fragment length.
#   map.qual: mapping quality to filter reads.
#   bowtie: boolean to indicate whether the aligner was Bowtie-like or not.
#   extr.only: boolean for doing extraction only.
#   strand.spec: string desc. for strand-specific coverage calculation.
# Return: list of coverage vectors, each vector represents a gene.

    strand.spec <- match.arg(strand.spec)
    # browser()
    g.start <- start(ranges(granges))  # vector of gene start sites.
    g.end <- end(ranges(granges))  # vector of gene end sites.
    gb.len <- g.end - g.start + 1
        
    # Filter gene ranges based on available chromosomes in bam file.
    # It seems scanBam will terminate if it querys a chromosome that does not
    # exist in bam file.
    inbam.mask <- as.vector(seqnames(granges)) %in% sn.inbam
    if(!any(inbam.mask)) {  # none matches.
        return(genZeroList(length(granges), gb.len))
    }

    # Scan bam file to retrieve short reads.
    sbw <- c('pos', 'qwidth', 'mapq', 'strand')  # scanBamWhat.
    sbp <- ScanBamParam(what=sbw, which=granges[inbam.mask], flag=scanBamFlag(
                        isUnmappedQuery=F, isNotPassingQualityControls=F, 
                        isDuplicate=F))
    sr.in.ranges <- tryCatch(scanBam(bam.file, param=sbp), error=function(e) e)
    if(class(sr.in.ranges)[1] == 'simpleError') {
        # This is not supposed to happen after those unmatched seqnames are
        # removed. I keep it for safety.
        return(genZeroList(length(granges), gb.len))
    }

    # Restore original order.
    sr.in.ranges <- sr.in.ranges[scanBamRevOrder(
                                 as.data.frame(granges[inbam.mask]), sbp)]

    # Calculate coverage for each gene range.
    scan.counter <- 0
    covg.allgenes <- vector('list', length(granges))  # coverage list to return.
    for(i in 1:length(granges)) {
        if(!inbam.mask[i] || extr.only) {
            covg.allgenes[[i]] <- Rle(0, gb.len[i])
            next
        }
        scan.counter <- scan.counter + 1
        # Retrieve short reads info as vectors.
        srg.pos <- sr.in.ranges[[scan.counter]]$pos
        srg.mapq <- sr.in.ranges[[scan.counter]]$mapq
        srg.strand <- sr.in.ranges[[scan.counter]]$strand
        srg.qwidth <- sr.in.ranges[[scan.counter]]$qwidth

        # Special handling for bowtie mapping.
        # browser()
        if(bowtie) {
            srg.mapq[is.na(srg.mapq)] <- 254
        }

        # Subset by mapping quality.
        q.mask <- which(srg.mapq >= map.qual)
        # q.mask <- which(srg.mapq >= map.qual)
        srg.pos <- srg.pos[q.mask]
        srg.strand <- srg.strand[q.mask]
        srg.qwidth <- srg.qwidth[q.mask]

        # Subset by strand info.
        # browser()
        if(strand.spec != 'both') {
            if(strand.spec == 'same') {
                s.mask <- which(srg.strand == as.character(v.strand[i]))
            } else {
                s.mask <- which(srg.strand != as.character(v.strand[i]))
            }
            srg.pos <- srg.pos[s.mask]
            srg.strand <- srg.strand[s.mask]
            srg.qwidth <- srg.qwidth[s.mask]
        }

        if(length(srg.pos) > 0) {
            # Adjust negative read positions for physical coverage.
            neg.idx <- srg.strand == '-'
            srg.pos[neg.idx] <- srg.pos[neg.idx] - fraglen + 
                                srg.qwidth[neg.idx]
            
            # Shift reads by subtracting start positions.
            srg.pos <- srg.pos - g.start[i] + 1
            
            # Calculate physical coverage on the whole genebody.
            covg.allgenes[[i]] <- coverage(IRanges(start=srg.pos, 
                                                   width=fraglen), 
                                           width=gb.len[i], method='sort')
        } else {
            covg.allgenes[[i]] <- Rle(0, gb.len[i])
        }
    }
    covg.allgenes
}


covBamExons <- function(granges.list, v.strand, bam.file, sn.inbam, fraglen, 
                        map.qual=20, bowtie=F, 
                        strand.spec=c('both', 'same', 'opposite')) {
# Extract coverage vectors from bam file for a list of transcripts of multiple
# exons.
# Args:
#   granges.list: list of GRanges objects representing genomic coordinates 
#                 to extract coverge from. Each GRanges object has multiple
#                 ranges representing exons and flanking regions. 
#   v.strand: vector of strand info.
#   bam.file: character string refers to the path of a bam file.
#   sn.inbam: vector of chromosome names in the bam file.
#   fraglen: fragment length.
#   map.qual: mapping quality to filter reads.
#   bowtie: boolean to indicate whether the aligner was Bowtie-like or not.
#   strand.spec: string desc. for strand-specific coverage calculation.
# Return: list of coverage vectors, each vector represents a transcript.

    strand.spec <- match.arg(strand.spec)
    # Obtain start and end sites for each gene.
    v.start <- vector('integer', length=length(granges.list))
    v.end <- vector('integer', length=length(granges.list))
    for(i in 1:length(granges.list)) {
        gr <- granges.list[[i]]
        v.start[i] <- start(ranges(gr))[1]
        v.end[i] <- tail(end(ranges(gr)), n=1)
    }
    dna.len <- v.end - v.start + 1

    # Obtain mRNA(including flanking) length for each gene.
    rna.len <- sapply(granges.list, function(g) {
            gr <- ranges(g)
            sum(end(gr) - start(gr) + 1)
        })

    # Filter transcripts whose chromosomes do not match bam file.
    trans.sn <- sapply(granges.list, function(t) {
            as.character(seqnames(t)[1])
        })
    inbam.mask <- trans.sn %in% sn.inbam
    if(!any(inbam.mask)) {  # none matches.
        return(genZeroList(length(granges.list), rna.len))
    }

    # Flatten the GRanges list for scanBam to process in batch.
    granges.combined <- do.call('c', granges.list[inbam.mask])
    sbw <- c('pos', 'qwidth', 'mapq', 'strand')  # scanBamWhat.
    sbp <- ScanBamParam(what=sbw, which=granges.combined, flag=scanBamFlag(
                        isUnmappedQuery=F, isNotPassingQualityControls=F, 
                        isDuplicate=F))

    # Scan bam file to retrieve short reads.
    sr.in.ranges <- tryCatch(scanBam(bam.file, param=sbp), error=function(e) e)
    if(class(sr.in.ranges)[1] == 'simpleError') {
        # This is not supposed to happen after those unmatched seqnames are
        # removed. I keep it for safty.
        return(genZeroList(length(granges.list), rna.len))
    }

    # Restore the original order.
    sr.in.ranges <- sr.in.ranges[
                        scanBamRevOrder(as.data.frame(granges.combined), sbp)
                    ]

    # Convert each list to data.fram for easy combine and retrieval.
    sr.in.ranges <- lapply(sr.in.ranges, as.data.frame)

    # Create a list of breaks for merging short reads.
    v.nranges <- sapply(granges.list[inbam.mask], length)
    brk.left <- 1
    brk.right <- -1
    v.brk.left <- vector('integer', length=length(granges.list[inbam.mask]))
    v.brk.right <- vector('integer', length=length(granges.list[inbam.mask]))
    for(i in 1:length(v.nranges)) {
        brk.right <- brk.left + v.nranges[i] - 1
        v.brk.left[i] <- brk.left
        v.brk.right[i] <- brk.right
        brk.left <- brk.right + 1
    }

    # Pool short reads that belong to the same transcript and 
    # calculate coverage.
    scan.counter <- 0
    covg.allgenes <- vector('list', length=length(granges.list))
    for(i in 1:length(granges.list)) {
        if(!inbam.mask[i]) {
            covg.allgenes[[i]] <- Rle(0, rna.len[i])
            next
        }
        scan.counter <- scan.counter + 1
        # Pool short reads from the same transcript.
        sr.pooled <- do.call('rbind', sr.in.ranges[v.brk.left[scan.counter]:
                                                   v.brk.right[scan.counter]])

        # Special handling for bowtie mapping.
        srg.mapq <- sr.pooled$mapq
        if(bowtie) {
            srg.mapq[is.na(srg.mapq)] <- 254
        }

        # Filter short reads by mapping quality.
        sr.pooled <- sr.pooled[which(srg.mapq >= map.qual), ]

        # Subset by strand info.
        if(strand.spec != 'both') {
            if(strand.spec == 'same') {
                s.mask <- which(sr.pooled$strand == as.character(v.strand[i]))
            } else {
                s.mask <- which(sr.pooled$strand != as.character(v.strand[i]))
            }
            sr.pooled <- sr.pooled[s.mask, ]
        }

        # Calculate coverage.
        if(nrow(sr.pooled) > 0) {
            # Adjust negative read positions for physical coverage.
            neg.idx <- sr.pooled$strand == '-'
            sr.pooled[neg.idx, ]$pos <- sr.pooled[neg.idx, ]$pos - fraglen + 
                                        sr.pooled[neg.idx, ]$qwidth
            
            # Shift reads by subtracting start positions.
            sr.pooled$pos <- sr.pooled$pos - v.start[i] + 1
            
            # Shift ranges by subtracting start positions.
            # BE careful with negative start positions! Need to adjust end
            # positions first(or the GRanges lib will emit errors if 
            # start > end). 
            # Negative start positions happen when flanking region exceeds 
            # the chromosomal start.
            gr <- ranges(granges.list[[i]])
            if(v.start[i] > 0) {
                start(gr) <- start(gr) - v.start[i] + 1
                end(gr) <- end(gr) - v.start[i] + 1
            } else {
                end(gr) <- end(gr) - v.start[i] + 1
                start(gr) <- start(gr) - v.start[i] + 1
            }
            
            # Calculate physical coverage on the whole genebody.
            covg <- coverage(IRanges(start=sr.pooled$pos, width=fraglen), 
                             width=dna.len[i], method='sort')
            
            # browser()
            # Concatenate all exon coverages.
            # covg.allgenes[[i]] <- seqselect(covg, gr)
            covg.allgenes[[i]] <- covg[gr]
        } else {
            covg.allgenes[[i]] <- Rle(0, rna.len[i])
        }
    }
    covg.allgenes
}

bamFileList <- function(ctg.tbl) {
# Determine the bam files involved in the configuration and whether it is a 
# bam file pair setup.
# Args:
#   ctg.tbl: coverage-genelist-title table.

    cov.uniq <- unique(ctg.tbl$cov)
    cov.list <- strsplit(cov.uniq, ':')
    v.nbam <- sapply(cov.list, length)
    v.bbp <- v.nbam == 2
    if(all(v.bbp)) {
        bbp <- T
    } else if(all(!v.bbp)) {
        bbp <- F
    } else {
        stop("No mix of bam and bam-pair allowed in configuration.\n")
    }

    list(bbp=bbp, bam.list=unique(unlist(cov.list)))
}

estiMapqStyle <- function(bam.file){
# Estimate the mapping quality style. Return TRUE if it is SAM standard.
# Sample 1000 mapped reads from bam file, and if the mapq of reads
# over half are NA, then return FALSE, because it is quite possible that
# the aligner using coding style as bowtie, 255 as highest score.
# Args:
#   bam.file: bam file to be sampled.

    sbw <- c('pos', 'qwidth', 'mapq', 'strand')
    sbp <- ScanBamParam(what=sbw, flag=scanBamFlag(
                        isUnmappedQuery=F, isDuplicate=F))
    samp <- BamSampler(bam.file, yieldSize=500)
    samp.reads <- scanBam(samp, param=sbp)[[1]]
    samp.len <- length(samp.reads[["mapq"]])
    mapq.255 <- sum(is.na(samp.reads[["mapq"]]))
    if(mapq.255/samp.len >= 0.5){
        return(FALSE)
    }else{
        return(TRUE)
    }
}

headerIndexBam <- function(bam.list) {
# Read bam header to determine mapping method.
# Index bam files if not done yet.
# Args:
#   ctg.tbl: coverage-genelist-title table.

    v.map.bowtie <- vector('logical', length=length(bam.list))
    for(i in 1:length(bam.list)) {
        bam.file <- bam.list[i]

        # Index bam file.
        if(!file.exists(paste(bam.file, ".bai", sep=""))) {
            indexBam(bam.file)
        }

        # Derive mapping program.
        header <- scanBamHeader(bam.file)
        map.prog <- try(strsplit(header[[1]]$text$'@PG'[[1]], ':')[[1]][2], 
                        silent=T)
        if(class(map.prog) != "try-error") {
            map.style <- grepl('tophat|bowtie|bedtools|star', map.prog, 
                                     ignore.case=T)
            if(map.style){
                v.map.bowtie[i] <- TRUE
                next
            }
            map.style <- grepl('bwa|casava', map.prog, 
                                     ignore.case=T)
            if(map.style){
                v.map.bowtie[i] <- FALSE
                next
            }
            if(estiMapqStyle(bam.file)){
                warning(sprintf("Aligner for: %s cannot be determined. Style of 
standard SAM mapping score will be used.", bam.file))
                v.map.bowtie[i] <- FALSE
            }else{
                warning(sprintf("Aligner for: %s cannot be determined. Style of 
Bowtie-like SAM mapping score will be used. Would you mind to tell us what 
aligner you are using?", bam.file))
                v.map.bowtie[i] <- TRUE
            }
        } else {
            cat("\n")
            if(estiMapqStyle(bam.file)){
                warning(sprintf("Aligner for: %s cannot be determined. Style of 
standard SAM mapping score will be used.", bam.file))
                v.map.bowtie[i] <- FALSE
            }else{
                warning(sprintf("Aligner for: %s cannot be determined. Style of 
Bowtie-like SAM mapping score will be used.", bam.file))
                v.map.bowtie[i] <- TRUE
            }
        }
    }
    names(v.map.bowtie) <- bam.list

    v.map.bowtie
}

libSizeBam <- function(bam.list) {
# Obtain library sizes by counting qualified bam records.
# Args:
#   ctg.tbl: coverage-genelist-title table.

    # Count only reads that are mapped, primary, passed quality control and 
    # un-duplicated.
    sbp <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F, isNotPrimaryRead=F, 
                        isNotPassingQualityControls=F, isDuplicate=F))
    v.lib.size <- vector('integer', length=length(bam.list))
    for(i in 1:length(bam.list)) {
        bfn <- bam.list[i]  # bam file name.
        cfn <- paste(basename(bfn), '.cnt', sep='')  # count file name.
        if(file.exists(cfn)) {
            v.lib.size[i] <- as.integer(readLines(cfn, n=1))
        } else {
            cnt.bam <- countBam(bfn, param=sbp)
            v.lib.size[i] <- cnt.bam$records
            writeLines(as.character(v.lib.size[i]), cfn)
        }
        names(v.lib.size)[i] <- bfn
    }

    v.lib.size
}

seqnamesBam <- function(bam.list) {
# Obtain chromosome names for each bam file. This list must be used to filter
# genomic regions before scanBam or it terminates immaturely.
#   ctg.tbl: coverage-genelist-title table.

    # browser()
    sn.list <- lapply(scanBamHeader(bam.list), function(h) {
            names(h$targets)
        })
    names(sn.list) <- bam.list

    sn.list
}

chrTag <- function(sn.inbam) {
# Check whether the chromosome name format in the bam file contains 'chr' or not.
# Args:
#   sn.inbam: seqnames in the bam file.

    n.chr <- length(grep('chr', sn.inbam))
    if(n.chr == 0) {
        chr.tag <- F
    } else if(n.chr == length(sn.inbam)) {
        chr.tag <- T
    } else {
        return("Inconsistent chromosome names in bam file. Check bam header.")
    }

    chr.tag
}

chunkIndex <- function(tot.gene, gcs) {
# Create chunk indices according to total number of genes and chunk size.
# Args:
#   tot.gene: total number of genes.
#   gcs: gene chunk size.

    nchk <- ceiling(tot.gene / gcs)  # number of chunks.
    chkidx.list <- vector('list', length=nchk)  # chunk indices list.
    chk.start <- 1  # chunk start.
    i.chk <- idiv(tot.gene, chunkSize=gcs)
    for(i in 1:nchk) {
        chk.size <- nextElem(i.chk)
        chkidx.list[[i]] <- c(chk.start, chk.start + chk.size - 1)
        chk.start <- chk.start + chk.size
    }

    chkidx.list
}

covMatrix <- function(chkidx.list, coord, rnaseq.gb, exonmodel, libsize, 
                      spit.dot=T, ...) {
# Function to generate a coverage matrix for all genes.
# Args:
#   chkidx.list: list of (start, end) indices for each chunk.
#   coord: dataframe of gene coordinates.
#   rnaseq.gb: boolean for RNA-seq genebody plot.
#   exonmodel: exon model data object.
#   libsize: total read count for this bam file.
#   spit.dot: boolean to control sptting '.' to consoles.
# Return: normalized coverage matrix for all genes, each row represents a gene.

    
    # Extract coverage and combine into a matrix.
    result.matrix <- foreach(chk=chkidx.list, .combine='rbind', 
                             .multicombine=T) %dopar% {
        if(spit.dot) {
            cat(".")
        }
        i <- chk[1]:chk[2]  # chunk: start -> end
        # If RNA-seq, retrieve exon ranges.
        if(rnaseq.gb) {
            exonranges.list <- unlist(exonmodel[coord[i, ]$tid])
        } else {
            exonranges.list <- NULL
        }
        doCov(coord[i, ], exonranges.list, ...)
    }

    # Floor negative values which are caused by spline.
    result.matrix[result.matrix < 0] <- 0
    result.matrix / libsize * 1e6  # normalize to RPM.


    ########### For debug #############
    # result.matrix <- matrix(0, nrow=nrow(coord), ncol=101)
    # for(c in 1:length(chkidx.list)) {
    #     chk <- chkidx.list[[c]]
    #     i <- chk[1]:chk[2]  # chunk: start -> end
    #     cat(".")
    #     # If RNA-seq, retrieve exon ranges.
    #     if(rnaseq.gb) {
    #         exonranges.list <- unlist(exonmodel[coord[i, ]$tid])
    #     } else {
    #         exonranges.list <- NULL
    #     }
    #     result.matrix[i, ] <- doCov(coord[i, ], exonranges.list, ...)
    # }
    # # Floor negative values which are caused by spline.
    # # browser()
    # result.matrix[result.matrix < 0] <- 0
    # result.matrix / libsize * 1e6  # normalize to RPM.
    ########### For debug #############
}


doCov <- function(coord.mat, exonranges.list, chr.tag, pint, reg2plot, 
                  flanksize, flankfactor, m.pts, f.pts, ...) {
# Extract coverage from bam file into a matrix. According to the parameter
# values, call corresponding functions.
# Args:
#   coord.mat: matrix of genomic coordinates to extract coverage.
#   exonranges.list: list of IRanges objects, each represents a group of exons.
#   pint: boolean of point interval.
#   reg2plot: string of region to plot.
#   flanksize: flanking region size.
#   flankfactor: flanking region factor.
#   m.pts: data points for middle interval.
#   f.pts: data points for flanking region.
# Return: matrix of interpolated coverage: each row represents a gene; each 
#         column represents a data point. Coverage from exons are concatenated.

    v.chrom <- coord.mat$chrom
    if(!chr.tag) {
        v.chrom <- sub('chr', '', v.chrom)
    }
    v.chrom <- as.factor(v.chrom)
    v.strand <- as.factor(coord.mat$strand)

    # Figure out interval region sizes and calculate flanking region sizes.
    if(!pint) {  # interval regions.
        if(!is.null(exonranges.list)) {  # RNA-seq
            if(flankfactor > 0) {
                v.fls <- sapply(exonranges.list, function(t) {
                            sum(end(t) - start(t) + 1) * flankfactor
                        })
            } else {
                v.fls <- rep(flanksize, length=nrow(coord.mat))
            }
            extrCovExons(v.chrom, exonranges.list, v.fls, v.strand, 
                         m.pts, f.pts, ...)
        } else {  # ChIP-seq with intervals.
            v.start <- coord.mat$start
            v.end <- coord.mat$end
            if(flankfactor > 0) {
                v.fls <- round((v.end - v.start + 1) * flankfactor)
            } else {
                v.fls <- rep(flanksize, length=nrow(coord.mat))
            }
            extrCov3Sec(v.chrom, v.start, v.end, v.fls, v.strand, 
                        m.pts, f.pts, ...)
        }
    } else {  # point center with flanking regions.
        v.midp <- vector('integer', length=nrow(coord.mat))
        for(r in 1:nrow(coord.mat)) {
            if(reg2plot == 'tss' && v.strand[r] == '+' || 
               reg2plot == 'tes' && v.strand[r] == '-') {
                # the left site is center.
                v.midp[r] <- coord.mat$start[r]
            } else {  # this also includes BED supplied point center.
                v.midp[r] <- coord.mat$end[r]
            }
        }
        # browser()
        extrCovMidp(v.chrom, v.midp, flanksize, v.strand, m.pts + f.pts*2, ...)
    }
}
