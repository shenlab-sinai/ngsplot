chromFormat <- function(crn, ...){
# Format chromosome names to our standard.
# Args:
#   crn: character vector of chromosome names.
# Returns: character vector of formatted chromosome names.

    crn <- sub('.fa$', '', crn)
    nochr.i <- grep('^chr', crn, invert=T)
    crn[nochr.i] <- paste('chr', crn[nochr.i], sep='')
    crn
}

FIScoreIntersect <- function(db.info, v.finfo){
# Further info score based on intersection with database FI columns.
# Used to select database files.
# Args: 
#   db.info: one record of annotation database table.
#   v.finfo: further info that is split into vector.
# Returns: a score for the intersection between database and further info.

    length(intersect(as.vector(db.info[c("FI.1", "FI.2", "FI.3")]), v.finfo))
}

SetupPlotCoord <- function(args.tbl, ctg.tbl, default.tbl, dbfile.tbl, progpath, 
                           genome, reg2plot, lgint, flanksize, bed.file, 
                           samprate) {
# Load genomic coordinates for plot based on the input arguments.
# Args:
#   args.tbl: input argument table
#   ctg.tbl: coverage-genelist-title table
#   default.tbl: default settings of databases and plot setting
#   dbfile.tbl: details of database files(.RData).
#   progpath: program path derived from NGSPLOT
#   genome: genome name, such as mm9, hg19.
#   reg2plot: tss, tes, genebody, etc.
#   lgint: boolean tag for large interval.
#   flanksize: flanking region size.
#   samprate: sampling rate

    # Subset using genome-region combination.
    key <- default.tbl$Genome == genome & default.tbl$Region == reg2plot
    if(sum(key) == 0) {
        stop("The combination of genome and region does not exist.\nYou may need to install the genome or the region does not exist yet.\n")
    }
    anno.parameters <- default.tbl[key, ]
    db.match.mask <- dbfile.tbl$Genome == genome & 
                     dbfile.tbl$Region == reg2plot
    anno.db.candidates <- dbfile.tbl[db.match.mask, ]
    
    # Database flavor.
    if('-D' %in% names(args.tbl)){  
        database <- as.character(args.tbl['-D'])
        database.allowed <- unique(anno.db.candidates$DB)
        stopifnot(database %in% database.allowed)
    }else{
        database <- as.character(anno.parameters$DefaultDB)
    }

    anno.db.candidates <- anno.db.candidates[anno.db.candidates$DB == database, ]

    prefix <- file.path(progpath, 'database', genome)

    # Further info to subset genomic regions.
    if('-F' %in% names(args.tbl)){
        finfo <- as.character(args.tbl['-F'])
    } else {
        finfo <- NULL
    }

    if(reg2plot!='bed'){
        Labs <- unlist(strsplit(as.character(anno.parameters$PointLab), ","))
        if(length(Labs)==1){
            pint <- TRUE
        }else{
            pint <- FALSE
        }

        if(!is.null(finfo)){  # use finfo to subset DB tables.
            v.finfo <- unlist(strsplit(finfo, ','))
            # Get the intersect of the v.finfo and FI columns of databases, 
            # and choose the best fit.
            fi.score <- apply(anno.db.candidates, 1, FIScoreIntersect, 
                              v.finfo=v.finfo)
            anno.db.candidates <- anno.db.candidates[fi.score == max(fi.score), ]
            # RNA-seq tag.
            if(reg2plot == 'genebody' && 'rnaseq' %in% v.finfo) {
                rnaseq.gb <- TRUE
            }else{
                rnaseq.gb <- FALSE
            }
        } else {
            rnaseq.gb <- FALSE
        }
        anno.db.candidates <- anno.db.candidates[order(anno.db.candidates$dbScore), ]
        f.load <- file.path(prefix, anno.db.candidates$"db.file"[1])

        if (!file.exists(f.load)){
            stop("The requested database file does not exist. You may have a corrupted database.\nConsider reinstalling the genome.\n")
            # cat("\nDownloading database:\n")
            # download.file(anno.db.candidates$"URL"[1], destfile=f.load, 
            #               method="curl")
            # stopifnot(file.exists(f.load))
        }
    } else if(reg2plot == 'bed') {
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
        # f.load <- paste(f.load)
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
        em.load <- file.path(prefix, paste(genome, database, 'exonmodel.RData', 
                                           sep='.'))
        load(em.load)  # load exon models into variable: exonmodel
        res$exonmodel <- exonmodel
    }

    res
}




