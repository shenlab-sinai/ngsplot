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

ConfigTbl <- function(args.tbl, fraglen){
# Create configuration table from "-C" argument.
# Args:
#   args.tbl: named vector of program arguments.
#   fraglen: fragment length.
# Returns: dataframe of configuration.

    covfile <- args.tbl['-C']

    if(length(grep('.txt$', covfile, ignore.case=T)) > 0) {  # config file.
        ctg.tbl <- read.table(covfile, sep="\t", colClasses='character', 
                              comment.char='#')
        if(ncol(ctg.tbl) < 3) {
            stop("Configuration file must contain at least 3 columns! Insufficient information provided.\n")
        }
        colnames(ctg.tbl)[1:3] <- c('cov', 'glist', 'title')
        if(ncol(ctg.tbl) >= 4) {
            colnames(ctg.tbl)[4] <- 'fraglen'
            fraglen.sp <- strsplit(ctg.tbl$fraglen, ":")
            if(!all(sapply(fraglen.sp, function(x) {
                        length(x) == 1 || length(x) == 2}))) {
                stop("Fragment length format must be X or X:Y; X and Y are integers.\n")
            }
            if(!all(as.integer(unlist(fraglen.sp)) > 0)) {
                stop("Fragment length must be positive integers! Check your configuration file.\n")
            }
        } else {
            ctg.tbl <- data.frame(ctg.tbl, fraglen=as.character(fraglen),
                                  stringsAsFactors=F)
        }
        if(ncol(ctg.tbl) >= 5) {
            colnames(ctg.tbl)[5] <- 'color'
            # Validate color specifications.
            col.validated <- col2rgb(ctg.tbl$color)
        } else {
            ctg.tbl <- data.frame(ctg.tbl, color=NA)
        }
        ctg.tbl

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
        data.frame(cov=covfile, glist=glist, title=title, 
                   fraglen=as.character(fraglen), color=NA, 
                   stringsAsFactors=F)
    }
}

# Global variables defined here.
.go.allowed <- c('total', 'max', 'prod', 'diff', 'hc', 'pca', 'none', 'hc2', 
                 'km')

setupVars <- function(args.tbl, anno.tbl, go.allowed){
# Setup variables from program arguments.
# Args:
#   args.tbl: named vector of program arguments.
#   anno.tbl: the database defaults table.
# Returns: list of variables.

    vl <- list()  # variable list to be exported.
    vl$genome <- args.tbl['-G']
    vl$reg2plot <- args.tbl['-R']
    vl$oname <- args.tbl['-O']

    #### Switch for Galaxy usage ####
    if('-Galaxy' %in% names(args.tbl)){
       stopifnot(as.integer(args.tbl['-Galaxy']) >= 0)
       vl$galaxy <- as.integer(args.tbl['-Galaxy'])
       vl$avgname <- args.tbl['-O2']
       vl$heatmapname <- args.tbl['-O3']
    }else{
       vl$galaxy <- as.integer(0)
    }
    
    #### Check region to plot ####
    region.allowed <- c(as.vector(unique(anno.tbl$Region)), "bed")
    if(!vl$reg2plot %in% region.allowed) {
        stop(paste(c("Unknown region specified. Must be one of:", 
                     region.allowed, "\n"), collapse=" "))
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
        flank.tbl <- c(2000, 2000, 2000, 500, 500, 1500, 1000, 1000)
        names(flank.tbl) <- c('tss','tes','genebody','exon','cgi', 'enhancer', 
                              'dhs','bed')
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

    #### Color scale string. ####
    if('-SC' %in% names(args.tbl)){ 
        vl$color.scale <- args.tbl['-SC']
        if(!vl$color.scale %in% c('local', 'region', 'global')) {
            scale.pair <- unlist(strsplit(vl$color.scale, ","))
            if(length(scale.pair) != 2 || is.na(as.numeric(scale.pair[1])) ||
               is.na(as.numeric(scale.pair[2]))) {
                stop("Color scale format error: must be local, region, global 
or a pair of numerics separated by ','\n")
            }
        }
    } else {
        vl$color.scale <- 'local'
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
        stopifnot(vl$go.algo %in% go.allowed)
    } else {
        vl$go.algo <- 'total'  # hierarchical clustering.
    }

    #### Algorithm for coverage vector normalization ####
    if('-AL' %in% names(args.tbl)) {
        vl$cov.algo <- args.tbl['-AL']
        al.allowed <- c('spline', 'bin')
        stopifnot(vl$cov.algo %in% al.allowed)
    } else {
        vl$cov.algo <- 'spline'
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
        stopifnot(vl$map.qual >= 0)
    } else {
        vl$map.qual <- 20
    }

    vl
}

replotVars <- function(args.tbl, existing.vl, bam.pair, go.allowed) {
# Setup replot variables.
# Args:
#   args.tbl: argument table.
#   existing.vl: existing variable list.
#   bam.pair: boolean for bam-pair use.

    vl <- list()  # variable list to be exported.
    # vl$iname <- args.tbl['-I']
    # vl$oname <- args.tbl['-O']

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

    #### Gene order algorithm ####
    if('-GO' %in% names(args.tbl)){ 
        vl$go.algo <- args.tbl['-GO']
        stopifnot(vl$go.algo %in% go.allowed)
    } 

    #### Reduce ratio. ####
    if('-RR' %in% names(args.tbl)){  
        vl$rr <- as.integer(args.tbl['-RR'])
        stopifnot(vl$rr >= 1)
    }

    # Variables that do not exist in ngs.plot.r or need a default value for
    # backward compatibility.

    #### Font size. ####
    if('-FS' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-FS']) > 0)
        vl$font.size <- as.integer(args.tbl['-FS'])
    } else {
        vl$font.size <- 12
    }

    #### Color scale string. ####
    if('-SC' %in% names(args.tbl)){ 
        vl$color.scale <- args.tbl['-SC']
        if(!vl$color.scale %in% c('local', 'region', 'global')) {
            scale.pair <- unlist(strsplit(vl$color.scale, ","))
            if(length(scale.pair) != 2 || is.na(as.numeric(scale.pair[1])) ||
               is.na(as.numeric(scale.pair[2]))) {
                stop("Color scale format error: must be local, region, global 
or a pair of numerics separated by ','\n")
            }
        }
    }

    #### Heatmap color. ####
    if('-CO' %in% names(args.tbl)) {
        vl$hm.color <- as.character(args.tbl['-CO'])
        v.colors <- unlist(strsplit(vl$hm.color, ":"))
        if(bam.pair && length(v.colors) != 2 || 
           !bam.pair && length(v.colors) != 1) {
            stop("Heatmap color specifications must correspond to bam-pair!\n")
        }
    } else {
        vl$hm.color <- "default"
    }

    # #### Remove zero tag. ####
    # if('-RZ' %in% names(args.tbl)){ 
    #     vl$rm.zero <- as.integer(args.tbl['-RZ'])
    #     stopifnot(vl$rm.zero == 0 || vl$rm.zero == 1)
    # } else if(!"rm.zero" %in% names(existing.vl)) {
    #     vl$rm.zero <- 1
    # }

    ##### Set cores number. ####
    if('-P' %in% names(args.tbl)){
        stopifnot(as.integer(args.tbl['-P']) >= 0)
        vl$cores.number <- as.integer(args.tbl['-P'])
    }else{
        vl$cores.number <- as.integer(0)
    }

    #### Misc. options for avg. profiles. ####
    vl$prof.misc <- list()
    if('-LEG' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-LEG']) >= 0)
        vl$prof.misc$legend <- as.integer(args.tbl['-LEG'])
    } else {
        vl$prof.misc$legend <- T
    }
    if('-BOX' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-BOX']) >= 0)
        vl$prof.misc$box <- as.integer(args.tbl['-BOX'])
    } else {
        vl$prof.misc$box <- T
    }
    if('-VLN' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-VLN']) >= 0)
        vl$prof.misc$vline <- as.integer(args.tbl['-VLN'])
    } else {
        vl$prof.misc$vline <- T
    }
    if('-XYL' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-XYL']) >= 0)
        vl$prof.misc$xylab <- as.integer(args.tbl['-XYL'])
    } else {
        vl$prof.misc$xylab <- T
    }
    if('-LWD' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-LWD']) > 0)
        vl$prof.misc$line.wd <- as.integer(args.tbl['-LWD'])
    } else {
        vl$prof.misc$line.wd <- 3
    }

    vl
}
############### End arguments configuration #####################
#################################################################
