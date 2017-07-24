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
    interflank.cov <- covBamExons(interflank.gr, v.strand, ...)

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
    cov.list <- covBamExons(granges, v.strand, ...)

    # Trim buffers from coverage vectors.
    cov.list <- lapply(cov.list, function(v) {
            window(v, start=bufsize + 1, width=length(v) - 2 * bufsize)
        })

    # Interpolate and reverse coverage vectors and assemble into a matrix.
    cov.mat <- matrix(nrow=length(cov.list), ncol=pts)
    for(i in 1:length(cov.list)) {
        if(is.null(cov.list[[i]])) {
            cov.mat[i, ] <- vector('integer', length=pts)
        } else {
            if(algo == "spline") {
                cov.mat[i, ] <- spline(1:length(cov.list[[i]]), cov.list[[i]], 
                                       n=pts)$y
            } else if(algo == "bin") {
                cov.mat[i, ] <- Bin(cov.list[[i]], n=pts)
            } else {
                # pass.
            }
            if(v.strand[i] == '-') {
                cov.mat[i, ] <- rev(cov.mat[i, ])
            }
        }
    }
    cov.mat
}

scanBamRevOrder <- function(org.gr, sbp) {
# ScanBamParam re-arranges the input genomic ranges. Use range info to 
# construct a string vector to find the order to reverse it.
    org.grnames <- with(org.gr, paste(seqnames, start, end, sep=':'))
    sbw.gr <- as.data.frame(bamWhich(sbp))  # scan-bam-ed
    if('space' %in% names(sbw.gr)) {
        sbw.grnames <- with(sbw.gr, paste(space, start, end, sep=':'))
    } else if('group_name' %in% names(sbw.gr)) {
        sbw.grnames <- with(sbw.gr, paste(group_name, start, end, sep=':'))
    } else {
        stop("Cannot locate chromosome names in extracted short reads. Report 
this problem using issue tracking or discussion forum.\n")
    }
    
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

covBamExons <- function(granges.dat, v.strand, bam.file, sn.inbam, fraglen, 
                        map.qual=20, bowtie=F, 
                        strand.spec=c('both', 'same', 'opposite')) {
# Extract coverage vectors from bam file for a list of transcripts of multiple
# exons.
# Args:
#   granges.dat: a GRanges object representing a set of genomic ranges or a 
#                list of GRanges objects each representing a set of exonic 
#                ranges. 
#   v.strand: vector of strand info.
#   bam.file: character string refers to the path of a bam file.
#   sn.inbam: vector of chromosome names in the bam file.
#   fraglen: fragment length.
#   map.qual: mapping quality to filter reads.
#   bowtie: boolean to indicate whether the aligner was Bowtie-like or not.
#   strand.spec: string desc. for strand-specific coverage calculation.
# Return: list of coverage vectors, each vector represents a transcript.

    strand.spec <- match.arg(strand.spec)

    if(class(granges.dat) == 'list') {
        # Construct a GRanges object representing DNA sequences.
        v.seqnames <- sapply(granges.dat, 
                             function(x) as.character(seqnames(x)[1]))
        v.start <- sapply(granges.dat, function(x) start(ranges(x))[1])
        v.end <- sapply(granges.dat, function(x) tail(end(ranges(x)), n=1))
        granges.dna <- GRanges(seqnames=v.seqnames, 
                               ranges=IRanges(start=v.start, end=v.end))
        # Obtain mRNA(including flanking) length for each gene.
        repr.lens <- sapply(granges.dat, function(g) {
            sum(end(g) - start(g) + 1)
        })
    } else {
        v.seqnames <- as.character(seqnames(granges.dat))
        v.start <- start(granges.dat)
        v.end <- end(granges.dat)
        granges.dna <- granges.dat
        repr.lens <- v.end - v.start + 1
        granges.dat <- vector('list', length(granges.dna)) # set null tags.
    }

    # Filter transcripts whose chromosomes do not match bam file.
    inbam.mask <- as.character(seqnames(granges.dna)) %in% sn.inbam
    if(!any(inbam.mask)) {  # none matches.
        return(genZeroList(length(granges.dna), repr.lens))
    }

    # scanBamWhat: the info that need to be extracted from a bam file.
    sbw <- c('pos', 'qwidth', 'mapq', 'strand', 'rname', 
             'mrnm', 'mpos', 'isize')
    sbp <- ScanBamParam(what=sbw, which=granges.dna[inbam.mask], 
                        flag=scanBamFlag(isUnmappedQuery=F, 
                                         isNotPassingQualityControls=F, 
                                         isDuplicate=F))

    # Scan bam file to retrieve short reads.
    sr.in.ranges <- tryCatch(scanBam(bam.file, param=sbp), error=function(e) e)
    if(class(sr.in.ranges)[1] == 'simpleError') {
        # This is not supposed to happen after those unmatched seqnames are
        # removed. I keep it for safty.
        return(genZeroList(length(granges.dna), repr.lens))
    }

    # Restore the original order.
    sr.in.ranges <- sr.in.ranges[scanBamRevOrder(
                                 as.data.frame(granges.dna[inbam.mask]), sbp)]

    CalcReadsCov <- function(srg, start, end, gr.rna, repr.len, strand) {
    # Calculate short read coverage for each gene/region.
    # Args:
    #   srg: extracted short reads in gene.
    #   start: start position of the DNA sequence.
    #   end: end position of the DNA sequence.
    #   gr.rna: GRanges object (multiple ranges) representing exon sequences.
    #           This can be NULL indicating the input ranges are DNAs.
    #   repr.len: DNA or mRNA sequence length.
    #   strand: transcript strand (+/-).
    # Returns: a coverage vector for the gene.

        # browser()
        # Special handling for bowtie mapping.
        if(bowtie) {
            srg <- within(srg, mapq[is.na(mapq)] <- 254)  # within!
        }
        # Filter short reads by mapping quality.
        all.mask <- srg$mapq >= map.qual

        # Subset by strand info.
        if(strand.spec != 'both') {
            if(strand.spec == 'same') {
                s.mask <- srg$strand == as.character(strand)
            } else {
                s.mask <- srg$strand != as.character(strand)
            }
            all.mask <- all.mask & s.mask
        }

        # If paired, filter reads that are not properly paired.
        paired <- all(with(srg, is.na(isize) | isize != 0))
        if(paired) {
            p.mask <- with(srg, rname == mrnm & xor(strand == '+', isize < 0))
            all.mask <- all.mask & p.mask
        }

        # Apply all the filters on short reads.
        srg <- lapply(srg, `[`, which(all.mask))

        # Calculate coverage.
        if(length(srg[[1]]) > 0) {
            if(paired) {
                cov.pos <- with(srg, ifelse(isize < 0, mpos, pos))
                cov.wd <- abs(srg$isize)
            } else {
                # Adjust negative read positions for physical coverage.
                cov.pos <- with(srg, ifelse(strand == '-', 
                                            pos - fraglen + qwidth, pos))
                cov.wd <- fraglen
            }
            # Shift reads by subtracting start positions.
            cov.pos <- cov.pos - start + 1
            # Calculate physical coverage on the whole genebody.
            covg <- coverage(IRanges(start=cov.pos, width=cov.wd), 
                             width=end - start + 1, method='sort')

            if(!is.null(gr.rna)) {  # RNA-seq.
                # Shift exonic ranges by subtracting start positions.
                # BE careful with negative start positions! Need to adjust end
                # positions first(or the GRanges lib will emit errors if 
                # start > end). 
                # Negative start positions happen when flanking region exceeds 
                # the chromosomal start.
                if(start > 0) {
                    start(gr.rna) <- start(gr.rna) - start + 1
                    end(gr.rna) <- end(gr.rna) - start + 1
                } else {
                    end(gr.rna) <- end(gr.rna) - start + 1
                    start(gr.rna) <- start(gr.rna) - start + 1
                }
                # Concatenate all exon coverages.
                covg[ranges(gr.rna)]
            } else {  # ChIP-seq.
                covg
            }
        } else {
            Rle(0, repr.len)
        }
    }

    covg.allgenes <- mapply(CalcReadsCov, srg=sr.in.ranges, 
                            start=v.start, end=v.end, gr.rna=granges.dat, 
                            repr.len=repr.lens, strand=v.strand, SIMPLIFY=F)
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
            map.style <- grepl('tophat|bowtie|bedtools|star|hisat', map.prog, 
                               ignore.case=T)
            if(map.style){
                v.map.bowtie[i] <- TRUE
                next
            }
            map.style <- grepl('bwa|casava|gem', map.prog, ignore.case=T)
            if(map.style) {
                v.map.bowtie[i] <- FALSE
                next
            }
        }
        warning(sprintf("Aligner for: %s cannot be determined. Style of 
standard SAM mapping score will be used. Would you mind submitting an issue 
report to us on Github? This will benefit people using the same aligner.", 
                        bam.file))
        v.map.bowtie[i] <- FALSE
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
    sbp <- ScanBamParam(flag=scanBamFlag(isUnmappedQuery=F, isSecondaryAlignment=F, 
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

    n.chr <- length(grep('^chr', sn.inbam))
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

covMatrix <- function(debug, chkidx.list, coord, rnaseq.gb, exonmodel, libsize, 
                      spit.dot=T, ...) {
# Function to generate a coverage matrix for all genes.
# Args:
#   debug: boolean tag for debugging.
#   chkidx.list: list of (start, end) indices for each chunk.
#   coord: dataframe of gene coordinates.
#   rnaseq.gb: boolean for RNA-seq genebody plot.
#   exonmodel: exon model data object.
#   libsize: total read count for this bam file.
#   spit.dot: boolean to control sptting '.' to consoles.
# Return: normalized coverage matrix for all genes, each row represents a gene.


    if(!debug) {
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

    } else {
        for(c in 1:length(chkidx.list)) {
            chk <- chkidx.list[[c]]
            i <- chk[1]:chk[2]  # chunk: start -> end
            cat(".")
            # If RNA-seq, retrieve exon ranges.
            if(rnaseq.gb) {
                exonranges.list <- unlist(exonmodel[coord[i, ]$tid])
            } else {
                exonranges.list <- NULL
            }
            cov <- doCov(coord[i, ], exonranges.list, ...)
            if(c == 1) {
                result.matrix <- matrix(0, nrow=nrow(coord), ncol=ncol(cov))
            }
            result.matrix[i, ] <- cov
        }
        # Floor negative values which are caused by spline.
        # browser()
        result.matrix[result.matrix < 0] <- 0
        result.matrix / libsize * 1e6  # normalize to RPM.

    }
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
    # if(!chr.tag) {
    #     v.chrom <- sub('chr', '', v.chrom)
    # }
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
