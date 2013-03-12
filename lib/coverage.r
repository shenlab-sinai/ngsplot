# Function to check if the range exceeds coverage vector boundaries.
checkBound <- function(start, end, range, chrlen){
    if(end + range > chrlen || start - range < 1)
        return(FALSE)   # out of boundary.
    else
        return(TRUE)
}

SplineRev3Sec <- function(cov.list, v.fls, pts.list, v.strand) {
# For 3 sections of continuous coverage, spline each according to specified 
# data points and return concatenated, interpolated curves.
# Args:
#   cov.list: a list of coverage vectors. Each vector represents a gene.
#   v.fls: vector of flanking region size.
#   pts.list: a list of three integers for data points.
#   v.strand: factor vector of strands.
# Return: matrix of interpolated coverage: each row represents a gene; each 
#         column represents a data point. Coverage from exons are concatenated.
# Notes: the names of cov.list and pts.list must be the same for index purpose.
# Suggested: use 'l', 'm', and 'r' for the names.

    # Create an empty coveage matrix first.
    tot.pts <- sum(unlist(pts.list))
    cov.mat <- matrix(nrow=length(cov.list), ncol=tot.pts)

    for(i in 1:length(cov.list)) {
        left.cov <- head(cov.list[[i]], n=v.fls[i])
        right.cov <- tail(cov.list[[i]], n=v.fls[i])
        mid.cov <- window(cov.list[[i]], start=v.fls[i] + 1, 
                          width=length(cov.list[[i]]) - 2*v.fls[i])
        left.cov <- spline(1:length(left.cov), left.cov, n=pts.list$l)$y
        right.cov <- spline(1:length(right.cov), right.cov, n=pts.list$r)$y
        mid.cov <- spline(1:length(mid.cov), mid.cov, n=pts.list$m)$y
        if(v.strand[i] == '+') {
            cov.mat[i, ] <- c(left.cov, mid.cov, right.cov)
        } else {
            cov.mat[i, ] <- rev(c(left.cov, mid.cov, right.cov))
        }
    }

    cov.mat
}

extrCov3Sec <- function(bam.file, sn.inbam, v.chrom, v.start, v.end, v.fls, 
                        bufsize, m.pts, f.pts, v.strand, fraglen, map.qual, 
                        bowtie) {
# Extract and interpolate coverage vectors from genomic regions with 3 sections.
# Args:
#   bam.file: character string refers to the path of a bam file.
#   v.chrom: factor vector of chromosome names.
#   v.start: integer vector of region start.
#   v.end: integer vector of region end.
#   v.fls: integer vector of flanking region size in bps.
#   m.pts: data points for middle interval.
#   f.pts: data points for flanking region.
#   v.strand: factor vector of gene strands.
#   fraglen: expected fragment length.
#   map.qual: map quality score cutoff.
# Return: matrix of interpolated coverage: each row represents a gene; each 
#         column represents a data point.

    interflank.gr <- GRanges(seqnames=v.chrom, ranges=IRanges(
                             start=v.start - v.fls - bufsize, 
                             end=v.end + v.fls + bufsize))
    interflank.cov <- covBam(bam.file, sn.inbam, interflank.gr, fraglen, 
                             map.qual, bowtie)

    # Trim buffers from coverage vectors.
    interflank.cov <- lapply(interflank.cov, function(v) {
            window(v, start=bufsize + 1, width=length(v) - 2 * bufsize)
        })

    # Interpolate and reverse coverage vectors.
    SplineRev3Sec(interflank.cov, v.fls, list(l=f.pts, m=m.pts, r=f.pts), 
                  v.strand)
}

# extrCovSec <- function(chrcov, start, end, flanking, pts, strand, weight){
# # This function should be deprecated.

#     cov <- as.vector(seqselect(chrcov, start - flanking, end + flanking))
#     cov.intp <- spline(1:length(cov), cov, n=pts)$y
#     if(weight){
#         cov.intp <- (length(cov) / pts) * cov.intp
#     }
#     if(strand == '+'){
#         cov.intp
#     }else{
#         rev(cov.intp)
#     }
# }


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

extrCovExons <- function(bam.file, sn.inbam, v.chrom, ranges.list, v.fls, 
                         bufsize, m.pts, f.pts, v.strand, fraglen, map.qual, 
                         bowtie) {
# Extract coverage vectors for transcripts with exon models.
# Args:
#   bam.file: character string refers to the path of a bam file.
#   v.chrom: factor vector of chromosome names.
#   ranges.list: list of IRanges objects for exon coordinates.
#   v.fls: integer vector of flanking region size in bps.
#   m.pts: data points for middle interval.
#   f.pts: data points for flanking region.
#   v.strand: factor vector of gene strands.
#   fraglen: expected fragment length.
#   map.qual: map quality score cutoff.
# Return: matrix of interpolated coverage: each row represents a gene; each 
#         column represents a data point. Coverage from exons are concatenated.

    # Construct ranges including exon and flanking regions.
    exonflank.list <- vector('list', length=length(ranges.list))
    for(i in 1:length(ranges.list)) {
        r.mod <- ranges.list[[i]]  # IRanges object to be modified.
        n <- length(r.mod)
        # Add flanking regions.
        start(r.mod)[1] <- start(r.mod)[1] - v.fls[i] - bufsize
        end(r.mod)[n] <- end(r.mod)[n] + v.fls[i] + bufsize
        exonflank.list[[i]] <- GRanges(seqnames=v.chrom[i], ranges=r.mod)
    }

    exonflank.cov <- covBamExons(bam.file, sn.inbam, exonflank.list, fraglen, 
                                 map.qual, bowtie)

    # Trim buffers from coverage vectors.
    exonflank.cov <- lapply(exonflank.cov, function(v) {
            window(v, start=bufsize + 1, width=length(v) - 2 * bufsize)
        })

    # Extract the left, exon and right parts into separate lists.
    # cov3.list <- sub3CovList(exonflank.cov, v.fls, v.fls)

    SplineRev3Sec(exonflank.cov, v.fls, list(l=f.pts, m=m.pts, r=f.pts), 
                  v.strand)
}


extrCovMidp <- function(bam.file, sn.inbam, v.chrom, v.midp, flanksize, bufsize, 
                        pts, v.strand, fraglen, map.qual, bowtie) {
# Extract coverage vectors with a middle point and symmetric flanking regions.
# Args:
#   bam.file: character string refers to the path of a bam file.
#   v.chrom: factor vector of chromosome names.
#   v.midp: integer vector of middle points.
#   flanksize: integer of flanking region size in bps.
#   pts: data points to spline into.
#   v.strand: factor vector of gene strands.
#   fraglen: expected fragment length.
#   map.qual: map quality score cutoff.
# Return: matrix of interpolated coverage: each row represents a gene; each 
#         column represents a data point.
    
    granges <- GRanges(seqnames=v.chrom, ranges=IRanges(
                       start=v.midp - flanksize - bufsize, 
                       end=v.midp + flanksize + bufsize))
    cov.list <- covBam(bam.file, sn.inbam, granges, fraglen, map.qual, bowtie)

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
            cov.mat[i, ] <- spline(1:length(cov.list[[i]]), cov.list[[i]], 
                                   n=pts)$y
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

covBam <- function(bam.file, sn.inbam, granges, fraglen, map.qual=20, bowtie=F) {
# Extract coverage vectors from bam file for a list of genes.
# Args:
#   bam.file: character string refers to the path of a bam file.
#   granges.list: list of GRanges objects representing genomic coordinates 
#                 to extract coverge from. In the case of RNA-seq, each GRanges
#                 object has multiple ranges representing exons and flanking 
#                 regions. In the case of ChIP-seq, each GRanges object has one
#                 range representing left boundary to right boundary.
# Return: list of coverage vectors, each vector represents a gene.

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
        # removed. I keep it for safty.
        return(genZeroList(length(granges), gb.len))
    }

    # Restore original order.
    sr.in.ranges <- sr.in.ranges[scanBamRevOrder(
                                 as.data.frame(granges[inbam.mask]), sbp)]

    # Calculate coverage for each gene range.
    scan.counter <- 0
    covg.allgenes <- vector('list', length(granges))  # coverage list to return.
    for(i in 1:length(granges)) {
        if(!inbam.mask[i]) {
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
        if(bowtie) {
            srg.mapq[is.na(srg.mapq)] <- 254
        }

        # Subset by mapping quality.
        q.mask <- srg.mapq >= map.qual
        srg.pos <- srg.pos[q.mask]
        srg.strand <- srg.strand[q.mask]
        srg.qwidth <- srg.qwidth[q.mask]

        if(length(srg.pos) > 0) {
            # Adjust negative read positions for physical coverage.
            neg.idx <- srg.strand == '-'
            srg.pos[neg.idx] <- srg.pos[neg.idx] - fraglen + srg.qwidth[neg.idx]
            
            # Shift reads by subtracting start positions.
            srg.pos <- srg.pos - g.start[i] + 1
            
            # Calculate physical coverage on the whole genebody.
            covg.allgenes[[i]] <- coverage(IRanges(start=srg.pos, width=fraglen), 
                                           width=gb.len[i], method='sort')
        } else {
            covg.allgenes[[i]] <- Rle(0, gb.len[i])
        }
    }
    covg.allgenes
}

covBamExons <- function(bam.file, sn.inbam, granges.list, fraglen, 
                        map.qual=20, bowtie=F) {
# Extract coverage vectors from bam file for a list of transcripts of multiple
# exons.
# Args:
#   bam.file: character string refers to the path of a bam file.
#   granges.list: list of GRanges objects representing genomic coordinates 
#                 to extract coverge from. Each GRanges object has multiple
#                 ranges representing exons and flanking regions. 
#   fraglen: fragment length.
#   map.qual: mapping quality to filter reads.
# Return: list of coverage vectors, each vector represents a transcript.

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
        if(bowtie) {
            sr.pooled$mapq[is.na(sr.pooled$mapq)] <- 254
        }

        # Filter short reads by mapping quality.
        sr.pooled <- sr.pooled[sr.pooled$mapq >= map.qual, ]

        # Calculate coverage.
        if(nrow(sr.pooled) > 0) {
            # Adjust negative read positions for physical coverage.
            neg.idx <- sr.pooled$strand == '-'
            sr.pooled[neg.idx, ]$pos <- sr.pooled[neg.idx, ]$pos - fraglen + 
                                        sr.pooled[neg.idx, ]$qwidth
            
            # Shift reads by subtracting start positions.
            sr.pooled$pos <- sr.pooled$pos - v.start[i] + 1
            
            # Shift ranges by subtracting start positions.
            gr <- ranges(granges.list[[i]])
            start(gr) <- start(gr) - v.start[i] + 1
            end(gr) <- end(gr) - v.start[i] + 1
            
            # Calculate physical coverage on the whole genebody.
            covg <- coverage(IRanges(start=sr.pooled$pos, width=fraglen), 
                             width=dna.len[i], method='sort')
            
            # Concatenate all exon coverages.
            covg.allgenes[[i]] <- seqselect(covg, gr)
        } else {
            covg.allgenes[[i]] <- Rle(0, rna.len[i])
        }
    }
    covg.allgenes
}

headerIndexBam <- function(ctg.tbl) {
# Read bam header to determine mapping method.
# Index bam files if not done yet.
# Args:
#   ctg.tbl: coverage-genelist-title table.

    cov.uniq <- unique(ctg.tbl$cov)
    v.map.bowtie <- vector('logical', length=length(cov.uniq))
    for(i in 1:length(cov.uniq)) {
        bam.file <- cov.uniq[i]

        # Index bam file.
        if(!file.exists(paste(bam.file, ".bai", sep=""))) {
            indexBam(bam.file)
        }

        # Derive mapping program.
        header <- scanBamHeader(bam.file)
        map.prog <- try(strsplit(header[[1]]$text$'@PG'[[1]], ':')[[1]][2])
        if(class(map.prog) != "try-error") {
            v.map.bowtie[i] <- ifelse(map.prog %in% c('Bowtie', 'TopHat'), T, F)
        } else {
            v.map.bowtie[i] <- F
        }
    }
    names(v.map.bowtie) <- cov.uniq

    v.map.bowtie
}

libSizeBam <- function(ctg.tbl) {
# Obtain library sizes by counting qualified bam records.
# Args:
#   ctg.tbl: coverage-genelist-title table.

    cov.uniq <- unique(ctg.tbl$cov)
    # Count only reads that are mapped, primary, passed quality control and 
    # un-duplicated.
    sbp <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F, isNotPrimaryRead=F, 
                        isNotPassingQualityControls=F, isDuplicate=F))
    v.lib.size <- vector('integer', length=length(cov.uniq))
    for(i in 1:length(cov.uniq)) {
        bfn <- cov.uniq[i]  # bam file name.
        cfn <- paste(bfn, '.cnt', sep='')  # count file name.
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

seqnamesBam <- function(ctg.tbl) {
# Obtain chromosome names for each bam file. This list must be used to filter
# genomic regions before scanBam or it terminates immaturely.
#   ctg.tbl: coverage-genelist-title table.

    cov.uniq <- unique(ctg.tbl$cov)
    sn.list <- lapply(scanBamHeader(cov.uniq), function(h) {
            bh.txt <- unlist(h$text)
            sn.txt <- bh.txt[grep('SN', bh.txt)]
            sub('SN:', '', sn.txt)
        })
    names(sn.list) <- basename(names(sn.list))

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

doCov <- function(coord.mat, chr.tag, reg2plot, pint, bam.file, sn.inbam, 
                  flanksize, flankfactor, bufsize, fraglen, map.qual, m.pts, 
                  f.pts, bowtie, exonranges.list=NULL) {
# Extract coverage from bam file into a matrix. According to the parameter
# values, call corresponding functions.
# Args:
#   coord.mat: matrix of genomic coordinates to extract coverage.
#   reg2plot: string of region to plot.
#   pint: boolean of point interval.
#   bam.file: character string refers to the path of a bam file.
#   flanksize: flanking region size.
#   flankfactor: flanking region factor.
#   fraglen: expected fragment length.
#   map.qual: map quality score cutoff.
#   m.pts: data points for middle interval.
#   f.pts: data points for flanking region.
#   bowtie: boolean of bowtie mapping method.
#   exonranges.list: list of IRanges objects, each represents a group of exons.
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
                            sum(end(t) - start(t) + 1)
                        })
            } else {
                v.fls <- rep(flanksize, length=nrow(coord.mat))
            }
            extrCovExons(bam.file, sn.inbam, v.chrom, exonranges.list, v.fls, 
                         bufsize, m.pts, f.pts, v.strand, fraglen, map.qual, 
                         bowtie)
        } else {  # ChIP-seq with intervals.
            v.start <- coord.mat$start
            v.end <- coord.mat$end
            if(flankfactor > 0) {
                v.fls <- round((v.end - v.start + 1) * flankfactor)
            } else {
                v.fls <- rep(flanksize, length=nrow(coord.mat))
            }
            extrCov3Sec(bam.file, sn.inbam, v.chrom, v.start, v.end, v.fls, 
                        bufsize, m.pts, f.pts, v.strand, fraglen, map.qual, 
                        bowtie)
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
        extrCovMidp(bam.file, sn.inbam, v.chrom, v.midp, flanksize, bufsize,
                    m.pts + f.pts*2, v.strand, fraglen, map.qual, bowtie)
    }
}
