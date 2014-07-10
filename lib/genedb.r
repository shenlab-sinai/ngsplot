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

SetupPlotCoord <- function(args.tbl, ctg.tbl, default.tbl, dbfile.tbl, 
                           progpath, genome, reg2plot, lgint, flanksize, 
                           samprate, galaxy) {
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
#   galaxy: boolean tag for galaxy installation.

    if(reg2plot!='bed'){
        # Subset using genome-region combination.
        key <- default.tbl$Genome == genome & default.tbl$Region == reg2plot
        if(sum(key) == 0) {
            stop("The combination of genome and region does not exist. You 
       may need to install the genome or the region does not exist yet.\n")
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

        if (file.exists(file.path(progpath, 'database', genome, reg2plot))){
            prefix <- file.path(progpath, 'database', genome, reg2plot)
        }else{
            prefix <- file.path(progpath, 'database', genome)
        }

        # Further info to subset genomic regions.
        if('-F' %in% names(args.tbl)){
            finfo <- as.character(args.tbl['-F'])
        } else {
            finfo <- NULL
        }

        Labs <- unlist(strsplit(as.character(anno.parameters$PointLab[1]), ","))
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
            stop("The requested database file does not exist. You may have a 
corrupted database. Consider reinstalling the genome.\n")
            # cat("\nDownloading database:\n")
            # download.file(anno.db.candidates$"URL"[1], destfile=f.load, 
            #               method="curl")
            # stopifnot(file.exists(f.load))
        }
        # Load genomic coordinates.
        cat("\nUsing database:\n")
        cat(paste(f.load, "\n", sep=""))
        load(f.load)  # load coordinates into variable: genome.coord

    } else if(reg2plot == 'bed') {
        rnaseq.gb <- FALSE
    }

    # Identify unique regions.
    reg.list <- as.factor(ctg.tbl$glist)
    uniq.reg <- levels(reg.list)

    # Determine plot coordinates for each unique region.
    coord.list <- vector('list', length(uniq.reg))
    names(coord.list) <- uniq.reg

    ReadBedCoord <- function(bed.file) {
    # Read a bed into memory as a dataframe.
    # Args:
    #   bed.file: path of the bed file to be read.
    # Returns: genome coordinates as a dataframe.

        if(!galaxy && 
            length(grep("\\.bed([0-9]+)?$", bed.file, ignore.case=T)) == 0) {
            warning(sprintf("File name: '%s' does not seem to a correct name 
for bed file.\n", bed.file))
        }
        bed.coord <- read.table(bed.file, sep="\t", quote="\"")
        if(ncol(bed.coord) <3){
            stop("A bed file must contain at least 3 columns! The format is: 
       chrom, start, end, gname, tid, strand. Columns 4-6 are optional\n")
        }
        genome.coord <- data.frame(chrom=chromFormat(bed.coord[, 1]), 
                                   start=bed.coord[, 2] + 1, 
                                   end=bed.coord[, 3], 
                                   gid=NA, gname='N', tid='N', strand='+', 
                                   byname.uniq=T, bygid.uniq=NA)
        # Perform sanity check for bed file.
        if(!all(genome.coord$start <= genome.coord$end)) {
            stop(sprintf("Sanity check for bed file: %s failed. Bed files are 
       0-based and right-side open.\n", bed.file))
        }

        # Deal with columns 4-6 (Optional).
        if(ncol(bed.coord) >=4){
            genome.coord$gname <- bed.coord[, 4]
        }
        if(ncol(bed.coord) >=5){
            genome.coord$tid <- bed.coord[, 5]
        }
        if(ncol(bed.coord) >=6){
            genome.coord$strand <- bed.coord[, 6]
        }

        genome.coord
    }

    # Sample a vector but preserve its original sequential order.
    sampleInSequence <- function(x, samprate) {
        x[ifelse(runif(length(x)) <= samprate, T, F)]
    }

    ValidateGeneSubset <- function(gid.match, tid.match, gname.match) {
    # To validate that there is no ambiguous hits in gene match.
    # Args:
    #   XXX.match: positional hits with respect to the gene database. 0 means 
    #              no hit.

        subset.matrix <- rbind(gid.match, tid.match, gname.match)
        g.valid <- sapply(subset.matrix, function(col.hits) {
            length(unique(col.hits[col.hits > 0])) <= 1
        })

        all(g.valid)
    }

    for(i in 1:length(uniq.reg)) {
        ur <- uniq.reg[i]
 
        if(reg2plot == "bed") {
            coord.list[[i]] <- ReadBedCoord(ur)
        } else if(ur == '-1') {  # use whole genome.
            coord.list[[i]] <- genome.coord[genome.coord$byname.uniq, ]
        } else {  # read gene list from text file.
            gene.list <- read.table(ur, as.is=T, comment.char='#')$V1
            gid.match <- match(gene.list, genome.coord$gid, nomatch=0)
            gname.match <- match(gene.list, genome.coord$gname, nomatch=0)
            tid.match <- match(gene.list, genome.coord$tid, nomatch=0)
            if(!ValidateGeneSubset(gid.match, tid.match, gname.match)) {
                stop("\nAmbiguous hits in gene database. This shall never 
       happen! Contact ngs.plot maintainers.\n")
            }
            subset.idx <- pmax(gid.match, tid.match, gname.match)
            # subset.idx <- gid.match + tid.match + gname.match
            subset.idx <- subset.idx[subset.idx != 0]
            if(length(subset.idx) == 0) {
                stop("\nGene subset size becomes zero. Are you using the 
       correct database?\n")
            }
            coord.list[[i]] <- genome.coord[subset.idx, ]
        }

        # Do sampling.
        if(samprate < 1) {
            samp.idx <- sampleInSequence(1:nrow(coord.list[[i]]), samprate)
            coord.list[[i]] <- coord.list[[i]][samp.idx, ]
        }
    }

    # If bed file, set boolean tag for point interval, i.e. interval=1bp.
    if(reg2plot == "bed") {
        pint <- all(sapply(coord.list, function(cd) all(cd$start == cd$end)))
        if(pint) {
            Labs <- c("Center")
        } else {
            Labs <- c("5'End", "3'End")
        }
    }

    # Set boolean tag for large interval display.
    if(!pint && is.na(lgint)) {
        if(median(sapply(coord.list, function(cd) {
                median(cd$end - cd$start + 1)
            })) > flanksize) {
            lgint <- 1
        } else {
            lgint <- 0
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




