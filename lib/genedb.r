# Function to remove ".fa" if needed, or to add "chr" if absent.
chromFormat <- function(crn, ...){
    crn <- sub('.fa$', '', crn)
    nochr.i <- grep('^chr', crn, invert=T)
    crn[nochr.i] <- paste('chr', crn[nochr.i], sep='')
    crn
}

filterFIScore <- function(db.info, finfos){
    score <- 3 - length(intersect(as.matrix(db.info[c("FI.1", "FI.2", "FI.3")]), as.matrix(finfos)))
    return(score)
}

SetupPlotCoord <- function(args.tbl, ctg.tbl, anno.tbl, anno.db.tbl, progpath, genome, reg2plot, 
                           lgint, flanksize, bed.file, samprate) {
# Load genomic coordinates for plot based on the input arguments.
# Args:
#   args.tbl: input argument table
#   ctg.tbl: coverage-genelist-title table
#   anno.tbl: default settings of databases and plot setting
#   anno.db.tbl: details of databases
#   progpath: program path derived from NGSPLOT
#   genome: genome name, such as mm9, hg19.
#   reg2plot: tss, tes, genebody, etc.
#   lgint: boolean tag for large interval.
#   flanksize: flanking region size.
#   samprate: sampling rate

    key <- paste(genome, reg2plot, sep='.')
    anno.parameters <- anno.tbl[key, ]
    # fi.allowed <- strsplit(as.character(anno.parameters$FurtherInfo), ',')[[1]]

    # Database flavor.
    if('-D' %in% names(args.tbl)){  
        database <- as.character(args.tbl['-D'])
        database.allowed <- as.character(unique(anno.db.tbl[which(anno.db.tbl$Genome==genome & anno.db.tbl$Region==reg2plot), "DB"]))
        stopifnot(database %in% database.allowed)
    }else{
        database <- anno.parameters$DefaultDB
    }

    prefix <- file.path(progpath, 'database', genome)

    # Further info to subset genomic regions.
    if('-F' %in% names(args.tbl)){
        finfo <- as.character(args.tbl['-F'])
    } else {
        finfo <- NULL
    }

    # Load coordinate dataset for the specified region. 
    # Use further info to subset regions.
    bed.tag <- F  # is the region a BED file?

    # if(reg2plot == 'tss' || reg2plot == 'tes') {
    #     f.load <- paste(prefix, 'genebody', sep='.')
    # } else if(reg2plot == 'genebody') {
    #     if(!is.null(finfo)){
    #         gb.allowed <- c('chipseq', 'rnaseq')
    #         stopifnot(finfo %in% gb.allowed)
    #     } else {
    #         finfo <- 'chipseq'            
    #     }
    #     f.load <- paste(prefix, 'genebody', sep='.')
    # } else if(reg2plot == 'exon') {
    #     if(!is.null(finfo)){
    #         exon.allowed <- c('canonical', 'variant', 'promoter', 'polyA', 
    #                           'altAcceptor', 'altDonor', 'altBoth')
    #         stopifnot(finfo %in% exon.allowed)
    #     } else {
    #         finfo <- 'canonical'            
    #     }
    #     f.load <- paste(prefix, 'exon', finfo, sep='.')
    # } else if(reg2plot == 'cgi') {
    #     if(!is.null(finfo)){
    #         cgi.allowed <- c("Genebody", "Genedesert", "OtherIntergenic", 
    #                          "Pericentromere", "Promoter1k", "Promoter3k", 
    #                          "ProximalPromoter")
    #         stopifnot(finfo %in% cgi.allowed)
    #     } else {
    #         finfo <- 'ProximalPromoter'
    #     }
    #     f.load <- paste(prefix, 'cgi', finfo, sep='.')
    # } else if(reg2plot == 'bed') {
    #     stopifnot(!is.null(bed.file))
    #     bed.coord <- read.table(bed.file, sep="\t")
    #     if(ncol(bed.coord) <3){
    #         stop('Input file must contain at least 3 columns!')
    #     }
    #     genome.coord <- data.frame(chrom=chromFormat(bed.coord[, 1]), 
    #                                start=bed.coord[, 2]+1, end=bed.coord[, 3], 
    #                                gid=NA, gname='N', tid='N', strand='+', 
    #                                byname.uniq=T, bygid.uniq=NA)
    #     if(ncol(bed.coord) >=4){
    #         genome.coord$gname <- bed.coord[, 4]
    #     }
    #     if(ncol(bed.coord) >=5){
    #         genome.coord$tid <- bed.coord[, 5]
    #     }
    #     if(ncol(bed.coord) >=6){
    #         genome.coord$strand <- bed.coord[, 6]
    #     }
    #     f.load <- NULL
    #     # bed.tag <- T
    # } else {
    #     # pass.
    # }

    if(reg2plot!='bed'){
        Labs <- strsplit(as.character(anno.parameters$PointLab), ",")[[1]]
        if(length(Labs)==1){
            pint <- TRUE
        }else{
            pint <- FALSE
        }
        anno.db.candidates <- anno.db.tbl[which(anno.db.tbl$Genome==genome & anno.db.tbl$Region==reg2plot), ]
        if(!is.null(finfo)){
            finfos <- strsplit(finfo, ',')[[1]]
            # get the intersect of the finfos and features of databases, and choose the best fit one
            anno.db.candidates$FIScore <- apply(anno.db.candidates, 1, filterFIScore, finfos=finfos)
            anno.db.candidates <- anno.db.candidates[which(anno.db.candidates$FIScore==anno.db.candidates[order(anno.db.candidates$FIScore), "FIScore"][1]), ]
            f.load <- file.path(prefix, anno.db.candidates[order(anno.db.candidates$dbScore), "db.file"][1])
            # RNA-seq tag.
            if(reg2plot == 'genebody' && ('rnaseq' %in% finfos)){
                rnaseq.gb <- TRUE
            }else{
                rnaseq.gb <- FALSE
            }
        }else{
            f.load <- file.path(prefix, anno.db.candidates[order(anno.db.candidates$dbScore), "db.file"][1])
            rnaseq.gb <- FALSE
        }
        if (!file.exists(f.load)){
            cat("\nDownloading database:\n")
            download.file(anno.db.candidates[order(anno.db.candidates$dbScore), "URL"][1], destfile=f.load, method="curl")
            stopifnot(file.exists(f.load))
        }
    }else if(reg2plot == 'bed') {
        stopifnot(!is.null(bed.file))
        bed.coord <- read.table(bed.file, sep="\t")
        if(ncol(bed.coord) <3){
            stop('Input file must contain at least 3 columns!')
        }
        genome.coord <- data.frame(chrom=chromFormat(bed.coord[, 1]), 
                                   start=bed.coord[, 2]+1, end=bed.coord[, 3], 
                                   gid=NA, gname='N', tid='N', strand='+', 
                                   byname.uniq=T, bygid.uniq=NA)
        if(ncol(bed.coord) >=4){
            genome.coord$gname <- bed.coord[, 4]
        }
        if(ncol(bed.coord) >=5){
            genome.coord$tid <- bed.coord[, 5]
        }
        if(ncol(bed.coord) >=6){
            genome.coord$strand <- bed.coord[, 6]
        }
        # Set tag for point interval, i.e. interval=1bp.
        if(all(genome.coord$end == genome.coord$start)) {
            pint <- TRUE
            Labs <- c("Center")
        } else{
            pint <- FALSE
            Labs <- c("Left", "Right")
        }
        f.load <- NULL
        rnaseq.gb <- FALSE
        # bed.tag <- T
    }

    # Load genomic coordinates.
    if(!is.null(f.load)) {
        cat("\nUsing database:\n")
        cat(paste(f.load, "\n", sep=""))
        f.load <- paste(f.load)
        load(f.load)  # load coordinates into variable: genome.coord
    }

    # Determine interval size automatically.
    if(!pint && is.na(lgint)) {
        if(median(genome.coord$end - genome.coord$start + 1) > flanksize) {
            lgint <- 1
        } else {
            lgint <- 0
        }
    }

    # Identify unique regions.
    reg.list <- as.factor(ctg.tbl$glist)
    uniq.reg <- levels(reg.list)

    # Determine plot coordinates for each unique region.
    coord.list <- vector('list', length(uniq.reg))
    names(coord.list) <- uniq.reg
    # Sample a vector but preserve its original sequential order.
    sampleInSequence <- function(x, samprate) {
        x[ifelse(runif(length(x)) <= samprate, T, F)]
    }
    for(i in 1:length(uniq.reg)) {
        ur <- uniq.reg[i]
 
        if(ur == '-1') {  # use whole genome.
            g.uniq <- which(genome.coord$byname.uniq)
            if(samprate < 1){
                coord.list[[i]] <- genome.coord[sampleInSequence(g.uniq, 
                                                                 samprate), ]
            } else {
                coord.list[[i]] <- genome.coord[g.uniq, ]
            }
        } else {  # read gene list from text file.
            gene.list <- read.table(ur, as.is=T, comment.char='#')$V1
            gid.match <- match(gene.list, genome.coord$gid, nomatch=0)
            gname.match <- match(gene.list, genome.coord$gname, nomatch=0)
            tid.match <- match(gene.list, genome.coord$tid, nomatch=0)
            subset.idx <- gid.match + tid.match + gname.match
            subset.idx <- subset.idx[subset.idx != 0]
            if(length(subset.idx) == 0) {
                stop("\nGene subset size becomes zero. Are you using the correct database?\n")
            }
            
            # subset.idx <- c(which(genome.coord$gname %in% gene.list & 
            #                       genome.coord$byname.uniq),
            #                 which(genome.coord$tid %in% gene.list))
            # # Test if all gid are NA. If database=refseq, this is true.
            # if(!all(is.na(genome.coord$gid))) {
            #     subset.idx <- c(subset.idx, 
            #                     which(genome.coord$gid %in% gene.list & 
            #                           genome.coord$bygid.uniq))
            # }
            # if(length(subset.idx) == 0) {
            #     stop("Gene subset size becomes zero. Are you using the correct database?\n")
            # }

            if(samprate < 1) {
                subset.idx <- sampleInSequence(subset.idx, samprate)
            }

            coord.list[[i]] <- genome.coord[subset.idx, ]
        }
    }

    res <- list(coord.list=coord.list, rnaseq.gb=rnaseq.gb, lgint=lgint,
                reg.list=reg.list, uniq.reg=uniq.reg, pint=pint, Labs=Labs)

    # Add exon models to the result if RNA-seq.
    if(rnaseq.gb) {
        em.load <- file.path(prefix, paste(genome, database, 'exonmodel.RData', sep='.'))
        load(em.load)  # load exon models into variable: exonmodel
        res$exonmodel <- exonmodel
    }

    res
}




