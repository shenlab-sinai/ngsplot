#!/usr/bin/env Rscript
#


# Deal with command line arguments.
cmd.help <- function(){
    cat("\nUsage: time.tabix.r -G genome -R region -C [cov|config]file\n")
    cat("                  [Options]\n")
    cat("\n## Mandatory parameters:\n")
    cat("  -G   Genome name, currently supported: mm9, hg19, rn4, sacCer3(ensembl only)\n")
    cat("  -R   Genomic regions to plot: tss, tes, genebody, exon, cgi or *.bed\n")
    cat("  -C   Tabix file\n")
    cat("## Optional parameters related to configuration file:\n")
    cat("  -E   Gene list to subset regions\n")
    cat("  -T   Image title\n")
    cat("## Important optional parameters:\n")
    cat("  -F   Further information can be provided to subset regions:\n")
    cat("         for genebody: chipseq(default), rnaseq.\n")
    cat("         for exon: canonical(default), variant, promoter, polyA,\n")
    cat("           altAcceptor, altDonor, altBoth.\n")
    cat("         for cgi: ProximalPromoter(default), Promoter1k, Promoter3k,\n")
    cat("           Genebody, Genedesert, OtherIntergenic, Pericentromere.\n")
    cat("  -D   Gene database: refseq(default), ensembl\n")
    cat("  -I   Shall interval be larger than flanking in plot?(0 or 1, default=automatic)\n")
    cat("  -L   Flanking region size\n")
    cat("  -N   Flanking region factor(will override flanking size)\n")
    cat("  -S   Randomly sample the regions for plot, must be:(0, 1]\n")
    cat("  -P   #CPUs to use. Set 0(default) for auto detection\n")
    # cat("  -PT  #data points to be used in plot(default=100)\n")
    cat("## Misc. parameters:\n")
    cat("  -GO  Gene order algorithm used in heatmaps: total(default), hc, max,\n")
    cat("         prod, diff, pca and none(according to gene list supplied)\n")
    cat("  -CS  Chunk size for loading genes in batch(default=100)\n")
    cat("  -FL  Fragment length used to calculate physical coverage(default=150)\n")
    cat("  -MQ  Mapping quality cutoff to filter reads(default=20)\n")
    cat("  -SE  Shall standard errors be plotted?(0 or 1)\n")
    cat("  -RB  The fraction of extreme values to be trimmed on both ends\n")
    cat("         default=0, 0.05 means 5% of extreme values will be trimmed\n")
    cat("  -FC  Flooding fraction:[0, 1), default=0.02\n")
    cat("  -FI  Forbid image output if set to 1(default=0)\n")
    cat("  -MW  Moving window width to smooth avg. profiles, must be integer\n")
    cat("         1=no(default); 3=slightly; 5=somewhat; 9=quite; 13=super.\n")
    cat("  -H   Opacity of shaded area, suggested value:[0, 0.5]\n")
    cat("         default=0, i.e. no shading, just curves\n")
    # cat("  -E   Calculate weighted coverage according to region size(experimental!)\n")
    cat("\n")
}


###########################################################################
#################### Deal with program input arguments ####################
args <- commandArgs(T)
# args <- unlist(strsplit('-G mm9 -R tss -C benchmark_morp_k9m2/k9_down_10M.bw -L 5000 -CS 1 -P 1 -S 0.1', ' '))

# Program environment variable.
progpath <- Sys.getenv('NGSPLOT')
if(progpath == "") {
    stop("Set environment variable NGSPLOT before run the program. \
          See README for details.\n")
}

# Input argument parser.
source(file.path(progpath, 'lib', 'parse.args.r'))
args.tbl <- parseArgs(args, c('-G', '-C', '-R'))
if(is.null(args.tbl)){
    cmd.help()
    stop('Error in parsing command line arguments. Stop.\n')
}

# Load required libraries.
if(!suppressMessages(require(ShortRead, warn.conflicts=F))) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(ShortRead)
    if(!suppressMessages(require(ShortRead, warn.conflicts=F))) {
        stop('Loading package ShortRead failed!')
    }
}
if(!suppressMessages(require(rtracklayer, warn.conflicts=F))) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(rtracklayer)
    if(!suppressMessages(require(rtracklayer, warn.conflicts=F))) {
        stop('Loading package rtracklayer failed!')
    }
}
if(!suppressMessages(require(BSgenome, warn.conflicts=F))) {
    source("http://bioconductor.org/biocLite.R")
    biocLite(BSgenome)
    if(!suppressMessages(require(BSgenome, warn.conflicts=F))) {
        stop('Loading package BSgenome failed!')
    }
}
if(!suppressMessages(require(doMC, warn.conflicts=F))) {
    install.packages('doMC')
    if(!suppressMessages(require(doMC, warn.conflicts=F))) {
        stop('Loading package doMC failed!')
    }
}
if(!suppressMessages(require(caTools, warn.conflicts=F))) {
    install.packages('caTools')
    if(!suppressMessages(require(caTools, warn.conflicts=F))) {
        stop('Loading package caTools failed!')
    }
}
if(!suppressMessages(require(utils, warn.conflicts=F))) {
    install.packages('utils')
    if(!suppressMessages(require(utils, warn.conflicts=F))) {
        stop('Loading package utils failed!')
    }
}

# Load tables of database: default.tbl, dbfile.tbl
default.tbl <- read.delim(file.path(progpath, 'database', 'default.tbl'))
dbfile.tbl <- read.delim(file.path(progpath, 'database', 'dbfile.tbl'))

# Setup variables from arguments.
argvar.list <- setupVars(args.tbl, default.tbl)
genome <- argvar.list$genome  # genome name, such as mm9, hg19, rn4.
reg2plot <- argvar.list$reg2plot  # tss, tes, genebody, bed...
bed.file <- argvar.list$bed.file  # BED file name if reg2plot=bed.
oname <- argvar.list$oname  # output file root name.
fi_tag <- argvar.list$fi_tag  # tag for forbidding image output
lgint <- argvar.list$lgint  # lgint: boolean tag for large interval
flanksize <- argvar.list$flanksize  # flanking region size
flankfactor <- argvar.list$flankfactor  # flanking region factor
samprate <- argvar.list$samprate  # sampling rate
shade.alp <- argvar.list$shade.alp  # shade area alpha
mw <- argvar.list$mw  # moving window width for smooth function.
cores.number <- argvar.list$cores.number  # #CPUs
se <- argvar.list$se  # se: tag for plotting stand errors
robust <- argvar.list$robust  # robust stat fraction
flood.frac <- argvar.list$flood.frac  # flooding fraction.
go.algo <- argvar.list$go.algo  # gene order algorithm used in heatmaps.
gcs <- argvar.list$gcs  # gcs: chunk size for grouping genes.
fraglen <- argvar.list$fraglen  # fragment length for physical coverage.
map.qual <- argvar.list$map.qual  # mapping quality cutoff.
bufsize <- argvar.list$bufsize  # buffer is to ensure smooth coverage at both 
                                # ends of the coverage vector.

# Configuration: coverage-genelist-title relationships.
ctg.tbl <- ConfigTbl(args.tbl, fraglen)

# Register doMC with CPU number.
if(cores.number == 0){
    registerDoMC()
} else {
    registerDoMC(cores.number)
}

# Setup plot-related coordinates and variables.
source(file.path(progpath, 'lib', 'genedb.r'))
plotvar.list <- SetupPlotCoord(args.tbl, ctg.tbl, default.tbl, dbfile.tbl, 
                               progpath, genome, reg2plot, lgint, flanksize, 
                               bed.file, samprate)
coord.list <- plotvar.list$coord.list  # list of coordinates for unique regions.
rnaseq.gb <- plotvar.list$rnaseq.gb  # tag for RNA-seq data.
lgint <- plotvar.list$lgint  # lgint: automatically determined if not specified.
reg.list <- plotvar.list$reg.list  # region list as in config file.
uniq.reg <- plotvar.list$uniq.reg  # unique region list.
pint <- plotvar.list$pint  # tag for point interval.
exonmodel <- plotvar.list$exonmodel  # exon ranges if rnaseq.gb=True.

# Setup data points for plot.
source(file.path(progpath, 'lib', 'plotlib.r'))
pts.list <- SetPtsSpline(pint, lgint)
pts <- pts.list$pts  # data points for avg. profile and standard errors.
m.pts <- pts.list$m.pts  # middle data points. For pint, m.pts=1.
f.pts <- pts.list$f.pts  # flanking region data points.

# Setup matrix for avg. profiles.
regcovMat <- CreatePlotMat(pts, ctg.tbl)
# Setup matrix for standard errors.
confiMat <- CreateConfiMat(se, pts, ctg.tbl)

# Genomic enrichment for all profiles in the config. Use this for heatmaps.
enrichList <- vector('list', nrow(ctg.tbl))

#######################################################################
# Here start to extract coverages for all genomic regions and calculate 
# data for plotting.

# Load coverage extraction lib.
source(file.path(progpath, 'lib', 'coverage.r'))

# Extract bam file names from configuration and determine if bam-pair is used.
bfl.res <- bamFileList(ctg.tbl)
bam.pair <- bfl.res$bbp  # boolean for bam-pair.
bam.list <- bfl.res$bam.list  # bam file list.

# Determine if bowtie is used to generate the bam files.
# Index bam files if not done yet.
# v.map.bowtie <- headerIndexBam(bam.list)

# Retrieve chromosome names for each bam file.
# sn.list <- seqnamesBam(bam.list)
seqnamesBW <- function(bam.list) {

    sn.list <- lapply(bam.list, function(t) {
            seqnames(seqinfo(BigWigFile(t)))
        })
    names(sn.list) <- bam.list

    sn.list
}
sn.list <- seqnamesBW(bam.list)

# Calculate library size from bam files for normalization.
# v.lib.size <- libSizeBam(bam.list)

############ Define some custome routines here ###############


covBW <- function(bw.file, sn.inbam, granges, fraglen, 
                  map.qual=20, bowtie=F) {

    # Wiggle is 1-based.
    g.start <- start(ranges(granges))  # vector of gene start sites.
    g.end <- end(ranges(granges))  # vector of gene end sites.
    gb.len <- g.end - g.start + 1
    chr.name <- seqnames(granges)
        
    # Filter gene ranges based on available chromosomes in bam file.
    # It seems scanBam will terminate if it querys a chromosome that does not
    # exist in bam file.
    inbam.mask <- as.vector(seqnames(granges)) %in% sn.inbam
    if(!any(inbam.mask)) {  # none matches.
        return(genZeroList(length(granges), gb.len))
    }

    # Retrieve coverage from bigWig.
    rec.in.ranges <- import(bw.file, which=granges[inbam.mask], 
                            asRangedData=F)

    # Trick: records for all query ranges are read into a GRanges object and 
    # genes cannot be distinguished from each other. So the gene coordinates 
    # have to be used to query the GRanges object in memory to retrieve the
    # records for a certain gene.

    # Calculate coverage for each gene range.
    lapply(1:length(granges), function(i) {
        v <- Rle(0, gb.len[i])
        if(!inbam.mask[i] || length(rec.in.ranges) == 0) {
            return(v)
        }
        # browser()
        rec.in.chrom <- rec.in.ranges[as.vector(seqnames(rec.in.ranges) == as.character(chr.name[i]))]
        if(length(rec.in.chrom) == 0) {
            return(v)
        }
        rec.in.gene <- rec.in.chrom[end(ranges(rec.in.chrom)) >= g.start[i] & 
                                    start(ranges(rec.in.chrom)) <= g.end[i]]
        if(length(rec.in.gene) == 0) {
            return(v)
        }
        r.starts <- start(ranges(rec.in.gene))
        r.ends <- end(ranges(rec.in.gene))
        s.starts <- r.starts - g.start[i] + 1
        s.ends <- r.ends - g.start[i] + 1
        s.starts[s.starts <= 1] <- 1
        s.ends[s.ends >= gb.len[i]] <- gb.len[i]

        # browser()
        for(j in 1:length(s.starts)) {
            v[s.starts[j]:s.ends[j]] <- score(rec.in.gene)[j]
        }

        v
    })
}

extrCov3SecBW <- function(bw.file, sn.inbam, v.chrom, v.start, v.end, v.fls, 
                        bufsize, m.pts, f.pts, v.strand, fraglen, map.qual, 
                        bowtie) {

    interflank.gr <- GRanges(seqnames=v.chrom, ranges=IRanges(
                             start=v.start - v.fls, 
                             end=v.end + v.fls))
    interflank.cov <- covBW(bw.file, sn.inbam, interflank.gr, fraglen, 
                             map.qual, bowtie)

    # Interpolate and reverse coverage vectors.
    SplineRev3Sec(interflank.cov, v.fls, list(l=f.pts, m=m.pts, r=f.pts), 
                  v.strand)
}

extrCovMidpBW <- function(bw.file, sn.inbam, v.chrom, v.midp, flanksize, 
                        bufsize, pts, v.strand, fraglen, map.qual, bowtie) {
    
    granges <- GRanges(seqnames=v.chrom, ranges=IRanges(
                       start=v.midp - flanksize, 
                       end=v.midp + flanksize))
    cov.list <- covBW(bw.file, sn.inbam, granges, fraglen, map.qual, bowtie)

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

doCovBW <- function(coord.mat, chr.tag, reg2plot, pint, bw.file, sn.inbam, 
                     flanksize, flankfactor, bufsize, fraglen, map.qual, m.pts, 
                     f.pts, bowtie, exonranges.list=NULL) {

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
            # extrCovExons(tbx.file, sn.inbam, v.chrom, exonranges.list, v.fls, 
            #              bufsize, m.pts, f.pts, v.strand, fraglen, map.qual, 
            #              bowtie)
        } else {  # ChIP-seq with intervals.
            v.start <- coord.mat$start
            v.end <- coord.mat$end
            if(flankfactor > 0) {
                v.fls <- round((v.end - v.start + 1) * flankfactor)
            } else {
                v.fls <- rep(flanksize, length=nrow(coord.mat))
            }
            extrCov3SecBW(bw.file, sn.inbam, v.chrom, v.start, v.end, v.fls, 
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
        extrCovMidpBW(bw.file, sn.inbam, v.chrom, v.midp, flanksize, bufsize,
                    m.pts + f.pts*2, v.strand, fraglen, map.qual, bowtie)
    }
}


covMatrixBW <- function(bw.file, libsize, sn.inbam, chr.tag, coord, 
                         chkidx.list, rnaseq.gb, exonmodel, reg2plot, pint, 
                         flanksize, flankfactor, bufsize, fraglen, map.qual, 
                         m.pts, f.pts, is.bowtie, spit.dot=T) {
    
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
        doCovBW(coord[i, ], chr.tag, reg2plot, pint, bw.file, sn.inbam, 
                 flanksize, flankfactor, bufsize, fraglen, map.qual, 
                 m.pts, f.pts, is.bowtie, exonranges.list)
    }

    # Floor negative values which are caused by spline.
    result.matrix[result.matrix < 0] <- 0

    # result.matrix / libsize * 1e6  # normalize to RPM.


    ########### For debug #############
    # pts <- m.pts + 2 * f.pts
    # result.matrix <- matrix(0, nrow=nrow(coord), ncol=pts)
    # for(c in 1:length(chkidx.list)) {
    #     chk <- chkidx.list[[c]]
    #     i <- chk[1]:chk[2]  # chunk: start -> end
    #     if(spit.dot) {
    #         cat(".")
    #     }
    #     # If RNA-seq, retrieve exon ranges.
    #     if(rnaseq.gb) {
    #         exonranges.list <- unlist(exonmodel[coord[i, ]$tid])
    #     } else {
    #         exonranges.list <- NULL
    #     }
    #     result.matrix[i, ] <- doCovBW(coord[i, ], chr.tag, reg2plot, 
    #           pint, bw.file, sn.inbam, flanksize, flankfactor, bufsize, 
    #           fraglen, map.qual, m.pts, f.pts, is.bowtie, 
    #           exonranges.list)
    # }
    # # Floor negative values which are caused by spline.
    # result.matrix[result.matrix < 0] <- 0

    # result.matrix / libsize * 1e6  # normalize to RPM.
    ########### For debug #############

    result.matrix
}




############ End of routines ############



# Process the config file row by row.
res.time <- system.time(for(r in 1:nrow(ctg.tbl)) {

    reg <- ctg.tbl$glist[r]  # retrieve gene list names.
    # Create coordinate chunk indices.
    chkidx.list <- chunkIndex(nrow(coord.list[[reg]]), gcs)

    # Do coverage for each bam file.
    bam.files <- unlist(strsplit(ctg.tbl$cov[r], ':'))
    # Obtain bam file basic info.
    # libsize <- v.lib.size[bam.files[1]]
    sn.inbam <- sn.list[[bam.files[1]]]
    chr.tag <- chrTag(sn.inbam)
    # is.bowtie <- v.map.bowtie[bam.files[1]]
    is.bowtie <- F
    if(class(chr.tag) == 'character') {
        stop(sprintf("Read %s error: %s", bam.files[1], chr.tag))
    }
    result.matrix <- covMatrixBW(bam.files[1], libsize, sn.inbam, chr.tag, 
                               coord.list[[reg]], chkidx.list, rnaseq.gb, 
                               exonmodel, reg2plot, pint, flanksize, 
                               flankfactor, bufsize, fraglen, map.qual, m.pts, 
                               f.pts, is.bowtie, spit.dot=F)
    if(bam.pair) {  # calculate background.
        pseudo.rpm <- 1e-9
        # libsize <- v.lib.size[bam.files[2]]
        sn.inbam <- sn.list[[bam.files[2]]]
        chr.tag <- chrTag(sn.inbam)
        # is.bowtie <- v.map.bowtie[bam.files[2]]
        is.bowtie <- F
        if(class(chr.tag) == 'character') {
            stop(sprintf("Read %s error: %s", bam.files[2], chr.tag))
        }
        bkg.matrix <- covMatrixBW(bam.files[2], libsize, sn.inbam, chr.tag, 
                                coord.list[[reg]], chkidx.list, rnaseq.gb, 
                                exonmodel, reg2plot, pint, flanksize, 
                                flankfactor, bufsize, fraglen, map.qual, m.pts, 
                                f.pts, is.bowtie, spit.dot=F)
        # browser()
        result.matrix <- log2((result.matrix + pseudo.rpm) / 
                              (bkg.matrix + pseudo.rpm))
    }

    # Calculate SEM if needed. Shut off SEM in single gene case.
    if(nrow(result.matrix) > 1 && se){
        confiMat[, r] <- apply(result.matrix, 2, function(x) CalcSem(x, robust))
    }

    # Book-keep this matrix for heatmap.
    enrichList[[r]] <- result.matrix

    # Return avg. profile.
    regcovMat[, r] <- apply(result.matrix, 2, function(x) mean(x, trim=robust, 
                                                               na.rm=T))
})

print(res.time)




