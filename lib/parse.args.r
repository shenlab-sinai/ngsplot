# Parse command line arguments and extract them into an associate array.
# Check if the required arguments are all satisfied.
parseArgs <- function(args, manditories){
    if(length(args) %% 2 == 1 || length(args) == 0){
        cat('Unpaired argument and value.\n')
        return(NULL)
    }
    n.i <- seq(1, length(args), by=2)
    v.i <- seq(2, length(args), by=2)
    args.name <- args[n.i]
    args.value <- args[v.i]

    # Check if required argument values are supplied.
    miss_tag <- F
    man.bool <- manditories %in% args.name
    if(!all(man.bool)){
        cat(paste('Missing argument: ', paste(manditories[!man.bool], 
                                              collapse=','), '.', sep='')
           )
        miss_tag <- T
    }
    if(miss_tag){
        res <- NULL
    }else{
        res <- args.value
        names(res) <- args.name
    }
    res
}

# Create configuration table from program arguments.
ConfigTbl <- function(args.tbl, ctg.tbl){

    covfile <- args.tbl['-C']

    if(length(grep('.txt$', covfile)) > 0) {  # config file.
        read.table(covfile, sep="\t", col.names=c('cov', 'glist', 'title'), 
                   colClasses='character', as.is=T, comment.char='#')
    } else {  # a single bam file.
        if('-E' %in% names(args.tbl)) {
            glist <- args.tbl['-E']
        } else {
            glist <- '-1'
        }
        if('-T' %in% names(args.tbl)) {
            title <- args.tbl['-T']
        } else {
            title <- 'Noname'
        }
        data.frame(cov=covfile, glist=glist, title=title, stringsAsFactors=F)
    }
}

# Setup misc. variables from program arguments for later use.
setupVars <- function(args.tbl, ctg.tbl){

    vl <- list()  # variable list to be exported.
    vl$genome <- args.tbl['-G']
    vl$reg2plot <- args.tbl['-R']
    vl$oname <- args.tbl['-O']

    prefix.regions <- c('tss','tes','genebody','exon','cgi')
    if(!vl$reg2plot %in% prefix.regions) {
        vl$bed.file <- vl$reg2plot  # save bed file name.
        vl$reg2plot <- 'bed'  # assume BED input.
    } else {
        vl$bed.file <- NULL
    }

    #### Image output forbidden tag. ####
    if('-FI' %in% names(args.tbl)){ 
        stopifnot(as.integer(args.tbl['-FI']) >= 0)
        vl$fi_tag <- as.integer(args.tbl['-FI'])
    }else{
        vl$fi_tag <- as.integer(0)
    }

    #### Create heatmap tag. ####
    # if('-HM' %in% names(args.tbl)){   
    #   stopifnot(as.integer(args.tbl['-HM']) >= 0)
    #   vl$hm_cols <- as.integer(args.tbl['-HM'])
    # }else{
    #   vl$hm_cols <- as.integer(100)
    # }

    #### Shall interval size be larger than flanking size? ####
    if('-I' %in% names(args.tbl)){  
        stopifnot(as.integer(args.tbl['-I']) >= 0)
        vl$lgint <- as.integer(args.tbl['-I'])
    }else{
        # int.tbl <- c(1, 0, 0, 1)
        # names(int.tbl) <- c('genebody', 'exon', 'cgi', 'bed')
        # vl$lgint <- int.tbl[vl$reg2plot]
        vl$lgint <- NA
    }

    #### Flanking region size. ####
    if('-L' %in% names(args.tbl)){  
        stopifnot(as.integer(args.tbl['-L']) >= 0)
        vl$flanksize <- args.tbl['-L']
    }else{
        flank.tbl <- c(2000, 2000, 2000, 500, 500, 1000)
        names(flank.tbl) <- c('tss','tes','genebody','exon','cgi','bed')
        vl$flanksize <- flank.tbl[vl$reg2plot]
    }
    vl$flanksize <- as.integer(vl$flanksize)

    #### Flanking size factor. ####
    if('-N' %in% names(args.tbl) && !('-L' %in% names(args.tbl))){  
        stopifnot(as.numeric(args.tbl['-N']) >= 0)
        vl$flankfactor <- as.numeric(args.tbl['-N'])
        # vl$flanksize <- floor(intsize*flankfactor)
    }else{
        vl$flankfactor <- 0.0
    }
    # if(vl$rnaseq.gb){ # RNA-seq plotting.
    #   vl$flanksize <- 0
    # }

    #### Random sampling rate. ####
    if('-S' %in% names(args.tbl)){  
        vl$samprate <- as.numeric(args.tbl['-S'])
        stopifnot(vl$samprate > 0 && vl$samprate <= 1)
    }else{
        vl$samprate <- 1.0
    }

    #### Shaded area alpha. ####
    if('-H' %in% names(args.tbl)){  
        vl$shade.alp <- as.numeric(args.tbl['-H'])
        stopifnot(vl$shade.alp >= 0 && vl$shade.alp < 1)
    }else{
        vl$shade.alp <- 0
    }

    #### Smooth function radius. ####
    if('-MW' %in% names(args.tbl)){  
        vl$mw <- as.integer(args.tbl['-MW'])
        stopifnot(vl$mw >= 1)
    }else{
        vl$mw <- 1
    }

    #### Weighted coverage. ####
    # if('-E' %in% names(args.tbl)){    
    #   stopifnot(as.integer(args.tbl['-E']) >= 0)
    #   vl$weight.genlen <- as.integer(args.tbl['-E'])
    # }else{
    #   vl$weight.genlen <- as.integer(0)
    # }

    ##### Set cores number. ####
    if('-P' %in% names(args.tbl)){
        stopifnot(as.integer(args.tbl['-P']) >= 0)
        vl$cores.number <- as.integer(args.tbl['-P'])
    }else{
        vl$cores.number <- as.integer(0)
    }

    #### Shall standard errors be plotted around average profiles? ####
    if('-SE' %in% names(args.tbl)){ 
        stopifnot(as.integer(args.tbl['-SE']) >= 0)
        vl$se <- as.integer(args.tbl['-SE'])
    } else{
        vl$se <- 1
    }

    #### Robust statistics. ####
    if('-RB' %in% names(args.tbl)){ 
        stopifnot(as.numeric(args.tbl['-RB']) >= 0)
        vl$robust <- as.numeric(args.tbl['-RB'])
    }else{
        vl$robust <- .0  # percentage to be trimmed on both ends.
    }

    #### Flooding fraction. ####
    if('-FC' %in% names(args.tbl)){ 
        vl$flood.frac <- as.numeric(args.tbl['-FC'])
        stopifnot(vl$flood.frac >= 0 && vl$flood.frac < 1)
    } else {
        vl$flood.frac <- .02  # <2% and >98% data will be illuminated equally.
    }

    #### Remove zero tag. ####
    if('-RZ' %in% names(args.tbl)){ 
        vl$rm.zero <- as.integer(args.tbl['-RZ'])
        stopifnot(vl$rm.zero == 0 || vl$rm.zero == 1)
    } else {
        vl$rm.zero <- 1  # remove all zero profiles in heatmaps.
    }

    #### Gene order algorithm ####
    if('-GO' %in% names(args.tbl)){ 
        vl$go.algo <- args.tbl['-GO']
        go.allowed <- c('total', 'max', 'prod', 'diff', 'hc', 'pca', 'none')
        stopifnot(vl$go.algo %in% go.allowed)
    } else {
        vl$go.algo <- 'total'  # hierarchical clustering.
    }

    #### Gene chunk size ####
    if('-CS' %in% names(args.tbl)) {
        vl$gcs <- as.integer(args.tbl['-CS'])
        stopifnot(vl$gcs > 0)
    } else {
        vl$gcs <- 100
    }

    #### Fragment length ####
    if('-FL' %in% names(args.tbl)) {
        vl$fraglen <- as.integer(args.tbl['-FL'])
        stopifnot(vl$fraglen > 0)
    } else {
        vl$fraglen <- 150
    }

    vl$bufsize <- vl$fraglen  # buffer added to both ends of the coverage vec.

    #### Mapping quality cutoff ####
    if('-MQ' %in% names(args.tbl)) {
        vl$map.qual <- as.integer(args.tbl['-MQ'])
        stopifnot(vl$map.qual > 0)
    } else {
        vl$map.qual <- 20
    }

    vl
}

replotVars <- function(args.tbl){

    vl <- list()  # variable list to be exported.
    vl$iname <- args.tbl['-I']
    vl$oname <- args.tbl['-O']

    #### Shaded area alpha. ####
    if('-H' %in% names(args.tbl)){  
        vl$shade.alp <- as.numeric(args.tbl['-H'])
        stopifnot(vl$shade.alp >= 0 && vl$shade.alp < 1)
    }

    #### Smooth function radius. ####
    if('-MW' %in% names(args.tbl)){  
        vl$mw <- as.integer(args.tbl['-MW'])
        stopifnot(vl$mw >= 1)
    }

    #### Shall standard errors be plotted around average profiles? ####
    if('-SE' %in% names(args.tbl)){ 
        stopifnot(as.integer(args.tbl['-SE']) >= 0)
        vl$se <- as.integer(args.tbl['-SE'])
    }

    #### Flooding fraction. ####
    if('-FC' %in% names(args.tbl)){ 
        vl$flood.frac <- as.numeric(args.tbl['-FC'])
        stopifnot(vl$flood.frac >= 0 && vl$flood.frac < 1)
    } 

    #### Heatmap color. ####
    if('-CO' %in% names(args.tbl)) {
        vl$hm.color <- as.character(args.tbl['-CO'])
    } else {
        vl$hm.color <- NULL
    }

    #### Remove zero tag. ####
    if('-RZ' %in% names(args.tbl)){ 
        vl$rm.zero <- as.integer(args.tbl['-RZ'])
        stopifnot(vl$rm.zero == 0 || vl$rm.zero == 1)
    }

    #### Gene order algorithm ####
    if('-GO' %in% names(args.tbl)){ 
        vl$go.algo <- args.tbl['-GO']
        go.allowed <- c('total', 'max', 'prod', 'diff', 'hc', 'pca', 'none')
        stopifnot(vl$go.algo %in% go.allowed)
    } 

    #### Reduce ratio. ####
    if('-RR' %in% names(args.tbl)){  
        vl$rr <- as.integer(args.tbl['-RR'])
        stopifnot(vl$rr >= 1)
    }

    ##### Set cores number. ####
    if('-P' %in% names(args.tbl)){
        stopifnot(as.integer(args.tbl['-P']) >= 0)
        vl$cores.number <- as.integer(args.tbl['-P'])
    }else{
        vl$cores.number <- as.integer(0)
    }

    vl
}
############### End arguments configuration #####################
#################################################################
