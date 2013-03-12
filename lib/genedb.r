# Function to remove ".fa" if needed, or to add "chr" if absent.
chromFormat <- function(crn, ...){
    crn <- sub('.fa$', '', crn)
    nochr.i <- grep('^chr', crn, invert=T)
    crn[nochr.i] <- paste('chr', crn[nochr.i], sep='')
    crn
}

SetupPlotCoord <- function(args.tbl, ctg.tbl, progpath, genome, reg2plot, 
                            samprate) {
# Load genomic coordinates for plot based on the input arguments.
# Args:
#   args.tbl: input argument table
#   ctg.tbl: coverage-genelist-title table
#   progpath: program path derived from NGSPLOT
#   genome: genome name, such as mm9, hg19.
#   reg2plot: tss, tes, genebody, etc.
#   samprate: sampling rate

    # Database flavor.
    if('-D' %in% names(args.tbl)){  
        database <- as.character(args.tbl['-D'])
        database.allowed <- c('refseq', 'ensembl')
        stopifnot(database %in% database.allowed)
    }else{
        database <- 'refseq'
    }

    prefix <- file.path(progpath, 'database', 
                        paste(genome, '.', database, sep=''))

    # Further info to subset genomic regions.
    if('-F' %in% names(args.tbl)){
        finfo <- as.character(args.tbl['-F'])
    } else {
        finfo <- NULL
    }

    # Load coordinate dataset for the specified region. 
    # Use further info to subset regions.
    bed.tag <- F  # is the region a BED file?

    if(reg2plot == 'tss' || reg2plot == 'tes') {
        f.load <- paste(prefix, 'genebody', sep='.')
    } else if(reg2plot == 'genebody') {
        if(!is.null(finfo)){
            gb.allowed <- c('chipseq', 'rnaseq')
            stopifnot(finfo %in% gb.allowed)
        } else {
            finfo <- 'chipseq'            
        }
        f.load <- paste(prefix, 'genebody', sep='.')
    } else if(reg2plot == 'exon') {
        if(!is.null(finfo)){
            exon.allowed <- c('canonical', 'variant', 'promoter', 'polyA', 
                              'altAcceptor', 'altDonor', 'altBoth')
            stopifnot(finfo %in% exon.allowed)
        } else {
            finfo <- 'canonical'            
        }
        f.load <- paste(prefix, 'exon', finfo, sep='.')
    } else if(reg2plot == 'cgi') {
        if(!is.null(finfo)){
            cgi.allowed <- c("Genebody", "Genedesert", "OtherIntergenic", 
                             "Pericentromere", "Promoter1k", "Promoter3k", 
                             "ProximalPromoter")
            stopifnot(finfo %in% cgi.allowed)
        } else {
            finfo <- 'ProximalPromoter'
        }
        f.load <- paste(prefix, 'cgi', finfo, sep='.')
    } else {  # assume reg2plot is a *.bed file.
        bed.coord <- read.table(reg2plot, sep="\t")
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
        f.load <- NULL
        bed.tag <- T
    }

    # Load genomic coordinates.
    if(!is.null(f.load)) {
        f.load <- paste(f.load, 'RData', sep='.')
        load(f.load)  # load coordinates into variable: genome.coord
    }

    # Identify unique regions.
    reg.list <- as.factor(ctg.tbl$glist)
    uniq.reg <- levels(reg.list)

    # Determine plot coordinates for each unique region.
    coord.list <- vector('list', length(uniq.reg))
    names(coord.list) <- uniq.reg
    for(i in 1:length(uniq.reg)) {
        ur <- uniq.reg[i]
 
        if(ur == '-1') {  # use whole genome.
            g.uniq <- which(genome.coord$byname.uniq)
            if(samprate < 1){
                samp.n <- round(samprate * length(g.uniq))
                coord.list[[i]] <- genome.coord[sample(g.uniq, samp.n), ]
            }else{
                coord.list[[i]] <- genome.coord[g.uniq, ]
            }
        } else {  # read gene list from text file.
            gene.list <- read.table(ur, as.is=T, comment.char='#')$V1
            subset.idx <- c(which(genome.coord$gname %in% gene.list & 
                                    genome.coord$byname.uniq),
                            which(genome.coord$tid %in% gene.list))
            # Test if all gid are NA. If database=refseq, this is true.
            if(!all(is.na(genome.coord$gid))) {
                subset.idx <- c(subset.idx, 
                                which(genome.coord$gid %in% gene.list & 
                                      genome.coord$bygid.uniq))
            }
            if(length(subset.idx) == 0) {
                stop("Gene subset size becomes zero. Are you using the correct database?\n")
            }
            if(samprate < 1) {
                samp.n <- round(samprate * length(subset.idx))
                subset.idx <- sample(subset.idx, samp.n)
            }
            coord.list[[i]] <- genome.coord[subset.idx, ]
        }
    }

    # Set tag for point interval, i.e. interval=1bp.
    if(reg2plot == 'tss' || reg2plot == 'tes' || bed.tag && 
        all(genome.coord$end == genome.coord$start)) {
        pint <- T
    } else{
        pint <- F
    }

    # RNA-seq tag.
    if(reg2plot == 'genebody' && finfo == 'rnaseq'){
        rnaseq.gb <- T
    }else{
        rnaseq.gb <- F
    }

    res <- list(coord.list=coord.list, rnaseq.gb=rnaseq.gb, 
        reg.list=reg.list, uniq.reg=uniq.reg, pint=pint)

    # Add exon models to the result if RNA-seq.
    if(rnaseq.gb) {
        em.load <- paste(prefix, 'exonmodel.RData', sep='.')
        load(em.load)  # load exon models into variable: exonmodel
        res$exonmodel <- exonmodel
    }

    res
}




