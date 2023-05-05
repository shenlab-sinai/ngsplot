# Parse command line arguments and extract them into an associate array.
# Check if the required arguments are all satisfied.
parseArgs <- function(args, manditories) {
    if(length(args) %% 2 == 1 || length(args) == 0) {
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

ConfigTbl <- function(args.tbl, fraglen) {
# Create configuration table from "-C" argument.
# Args:
#   args.tbl: named vector of program arguments.
#   fraglen: fragment length.
# Returns: dataframe of configuration.

    covfile <- args.tbl['-C']

    suppressWarnings(
        ctg.tbl <- tryCatch(
            read.table(covfile, colClasses='character', comment.char='#'), 
            error=function(e) {
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
                           fraglen=as.character(fraglen), 
                           color=NA, stringsAsFactors=F)
            }
        )
    )

    # Read a config file.
    if(ncol(ctg.tbl) < 3) {
        stop("Configuration file must contain at least 3 columns! 
Insufficient information provided.\n")
    }
    colnames(ctg.tbl)[1:3] <- c('cov', 'glist', 'title')
    if(ncol(ctg.tbl) >= 4) {
        colnames(ctg.tbl)[4] <- 'fraglen'
        fraglen.sp <- strsplit(ctg.tbl$fraglen, ":")
        if(!all(sapply(fraglen.sp, function(x) {
                    length(x) == 1 || length(x) == 2}))) {
            stop("Fragment length format must be X or X:Y; X and Y are 
integers.\n")
        }
        if(!all(as.integer(unlist(fraglen.sp)) > 0)) {
            stop("Fragment length must be positive integers! Check your 
configuration file.\n")
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
}

CheckRegionAllowed <- function(reg2plot, anno.tbl) {
# Check if region to plot is an allowed value.

    region.allowed <- c(as.vector(unique(anno.tbl$Region)), "bed")
    if(!reg2plot %in% region.allowed) {
        stop(paste(c("Unknown region specified. Must be one of:", 
                     region.allowed, "\n"), collapse=" "))
    }
}

CoverageVars <- function(args.tbl, reg2plot) {
# Setup variables from program arguments.
# Args:
#   args.tbl: named vector of program arguments.
#   reg2plot: string describing region to plot.
# Returns: list of variables.

    vl <- list()  # vl: configured variable list

    #### Switch for debug ####
    if('-Debug' %in% names(args.tbl)) {
       stopifnot(as.integer(args.tbl['-Debug']) >= 0)
       vl$debug <- as.integer(args.tbl['-Debug'])
    } else {
       vl$debug <- as.integer(0)
    }

    #### Switch for Galaxy usage ####
    if('-Galaxy' %in% names(args.tbl)) {
       stopifnot(as.integer(args.tbl['-Galaxy']) >= 0)
       vl$galaxy <- as.integer(args.tbl['-Galaxy'])
       vl$avgname <- args.tbl['-O2']
       vl$heatmapname <- args.tbl['-O3']
    } else {
       vl$galaxy <- as.integer(0)
    }
    
    ######## Coverage-generation parameters ########
    #### Flanking region size. ####
    if('-L' %in% names(args.tbl)){  
        vl$flanksize <- as.integer(args.tbl['-L'])
        stopifnot(vl$flanksize >= 0)
    } else {
        flank.tbl <- setNames(
            c(2000, 2000, 2000, 500, 500, 1500, 1000, 1000),
            c('tss','tes','genebody','exon','cgi', 'enhancer', 'dhs','bed')) 
        vl$flanksize <- as.integer(flank.tbl[reg2plot])
    }

    #### Flanking size factor. ####
    if('-N' %in% names(args.tbl) && !('-L' %in% names(args.tbl))) {
        stopifnot(as.numeric(args.tbl['-N']) >= 0)
        vl$flankfactor <- as.numeric(args.tbl['-N'])
    } else {
        vl$flankfactor <- 0.0
    }

    #### Robust statistics. ####
    if('-RB' %in% names(args.tbl)){ 
        stopifnot(as.numeric(args.tbl['-RB']) >= 0)
        vl$robust <- as.numeric(args.tbl['-RB'])
    }else{
        vl$robust <- .0  # percentage to be trimmed on both ends.
    }

    #### Random sampling rate. ####
    if('-S' %in% names(args.tbl)){  
        vl$samprate <- as.numeric(args.tbl['-S'])
        stopifnot(vl$samprate > 0 && vl$samprate <= 1)
    }else{
        vl$samprate <- 1.0
    }

    ##### Set cores number. ####
    if('-P' %in% names(args.tbl)){
        stopifnot(as.integer(args.tbl['-P']) >= 0)
        vl$cores.number <- as.integer(args.tbl['-P'])
    }else{
        vl$cores.number <- as.integer(0)
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

    #### Mapping quality cutoff ####
    if('-MQ' %in% names(args.tbl)) {
        vl$map.qual <- as.integer(args.tbl['-MQ'])
        stopifnot(vl$map.qual >= 0)
    } else {
        vl$map.qual <- 20
    }

    #### Fragment length ####
    if('-FL' %in% names(args.tbl)) {
        vl$fraglen <- as.integer(args.tbl['-FL'])
        stopifnot(vl$fraglen > 0)
    } else {
        vl$fraglen <- 150
    }
    vl$bufsize <- vl$fraglen  # buffer added to both ends of the coverage vec.

    #### Strand-specific coverage ####
    if('-SS' %in% names(args.tbl)) {
        spec.allowed <- c('both', 'same', 'opposite')
        stopifnot(args.tbl['-SS'] %in% spec.allowed)
        vl$strand.spec <- args.tbl['-SS']
    } else {
        vl$strand.spec <- 'both'
    }

    #### Shall interval size be larger than flanking size? ####
    if('-IN' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-IN']) >= 0)
        vl$inttag <- as.integer(args.tbl['-IN'])
    } else {
        vl$inttag <- NA
    }

    #### Image output forbidden tag. ####
    if('-FI' %in% names(args.tbl)){ 
        stopifnot(as.integer(args.tbl['-FI']) >= 0)
        vl$fi_tag <- as.integer(args.tbl['-FI'])
    }else{
        vl$fi_tag <- as.integer(0)
    }

    vl
}

PlotVars <- function(args.tbl, existing.vl=vector('character'), 
                     prof.misc=list(), low.count=NULL, go.paras=list()) {
# Setup replot variables.
# Args:
#   args.tbl: argument table.
#   existing.vl: existing variable name character list.
#   prof.misc: misc. avg prof variable list.
#   go.paras: gene ordering parameters.
# Returns: list of updated variables.

    ## Plotting-related parameters:
    updated.vl <- list()
    ### Misc. parameters:
    #### Font size. ####
    if('-FS' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-FS']) > 0)
        updated.vl$font.size <- as.integer(args.tbl['-FS'])
    } else if(!'font.size' %in% existing.vl) {
        updated.vl$font.size <- 12
    }

    ### Avg. profiles parameters:
    #### Plot width. ####
    if('-WD' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-WD']) > 0)
        updated.vl$plot.width <- as.integer(args.tbl['-WD'])
    } else if(!'plot.width' %in% existing.vl) {
        updated.vl$plot.width <- 8
    }

    #### Plot height. ####
    if('-HG' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-HG']) > 0)
        updated.vl$plot.height <- as.integer(args.tbl['-HG'])
    } else if(!'plot.height' %in% existing.vl) {
        updated.vl$plot.height <- 7
    }

    #### Shall standard errors be plotted around average profiles? ####
    if('-SE' %in% names(args.tbl)) { 
        stopifnot(as.integer(args.tbl['-SE']) >= 0)
        updated.vl$se <- as.integer(args.tbl['-SE'])
    } else if(!'se' %in% existing.vl) {
        updated.vl$se <- 1
    }

    #### Smooth function radius. ####
    if('-MW' %in% names(args.tbl)) {  
        stopifnot(as.integer(args.tbl['-MW']) >= 1)
        updated.vl$mw <- as.integer(args.tbl['-MW'])
    } else if(!'mw' %in% existing.vl) {
        updated.vl$mw <- 1
    }

    #### Shaded area alpha. ####
    if('-H' %in% names(args.tbl)) {
        stopifnot(as.numeric(args.tbl['-H']) >= 0 && 
                  as.numeric(args.tbl['-H']) < 1)
        updated.vl$shade.alp <- as.numeric(args.tbl['-H'])
    } else if(!'shade.alp' %in% existing.vl) {
        updated.vl$shade.alp <- 0
    }

    #### Misc. options for avg. profiles. ####
    updated.vl$prof.misc <- prof.misc
    if('-YAS' %in% names(args.tbl)) {
        ystr <- args.tbl['-YAS']
        if(ystr != 'auto') {
            yp <- unlist(strsplit(ystr, ','))
            if(length(yp) == 2) {
                y.min <- as.numeric(yp[1])
                y.max <- as.numeric(yp[2])
                stopifnot(y.min < y.max)
            } else {
                stop("-YAS must be 'auto' or a pair of numerics separated 
by ','\n")
            }
            updated.vl$prof.misc$yscale <- c(y.min, y.max)
        } else {
            updated.vl$prof.misc$yscale <- 'auto'
        }
    } else if(!'yscale' %in% names(prof.misc)) {
        updated.vl$prof.misc$yscale <- 'auto'
    }
    if('-LEG' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-LEG']) >= 0)
        updated.vl$prof.misc$legend <- as.integer(args.tbl['-LEG'])
    } else if(!'legend' %in% names(prof.misc)) {
        updated.vl$prof.misc$legend <- T
    }
    if('-BOX' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-BOX']) >= 0)
        updated.vl$prof.misc$box <- as.integer(args.tbl['-BOX'])
    } else if(!'box' %in% names(prof.misc)) {
        updated.vl$prof.misc$box <- T
    }
    if('-VLN' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-VLN']) >= 0)
        updated.vl$prof.misc$vline <- as.integer(args.tbl['-VLN'])
    } else if(!'vline' %in% names(prof.misc)) {
        updated.vl$prof.misc$vline <- T
    }
    if('-XYL' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-XYL']) >= 0)
        updated.vl$prof.misc$xylab <- as.integer(args.tbl['-XYL'])
    } else if(!'xylab' %in% names(prof.misc)) {
        updated.vl$prof.misc$xylab <- T
    }
    if('-LWD' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-LWD']) > 0)
        updated.vl$prof.misc$line.wd <- as.integer(args.tbl['-LWD'])
    } else if(!'line.wd' %in% names(prof.misc)) {
        updated.vl$prof.misc$line.wd <- 3
    }

    ### Heatmap parameters:
    #### Gene order algorithm ####
    if('-GO' %in% names(args.tbl)){ 
        go.allowed <- c('total', 'max', 'prod', 'diff', 'hc', 'none', 'km')
        stopifnot(args.tbl['-GO'] %in% go.allowed)
        updated.vl$go.algo <- args.tbl['-GO']
    } else if(!'go.algo' %in% existing.vl){
        updated.vl$go.algo <- 'total'  # hierarchical clustering.
    }

    #### Reduce ratio to control a heatmap height ####
    if('-RR' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-RR']) > 0)
        updated.vl$rr <- as.integer(args.tbl['-RR'])
    } else if(!'rr' %in% existing.vl) {
        updated.vl$rr <- 30
    }

    #### Color scale string. ####
    if('-SC' %in% names(args.tbl)) {
        if(!args.tbl['-SC'] %in% c('local', 'region', 'global')) {
            scale.pair <- unlist(strsplit(args.tbl['-SC'], ","))
            if(length(scale.pair) != 2 || is.na(as.numeric(scale.pair[1])) ||
               is.na(as.numeric(scale.pair[2]))) {
                stop("Color scale format error: must be local, region, global 
or a pair of numerics separated by ','\n")
            }
        }
        updated.vl$color.scale <- args.tbl['-SC']
    } else if(!'color.scale' %in% existing.vl) {
        updated.vl$color.scale <- 'local'
    }

    #### Flooding fraction. ####
    if('-FC' %in% names(args.tbl)) { 
        stopifnot(as.numeric(args.tbl['-FC']) >= 0 && 
                  as.numeric(args.tbl['-FC']) < 1)
        updated.vl$flood.frac <- as.numeric(args.tbl['-FC'])
    } else if(!'flood.frac' %in% existing.vl) {
        updated.vl$flood.frac <- .02
    }

    #### Heatmap color. ####
    if('-CO' %in% names(args.tbl)) {
        updated.vl$hm.color <- as.character(args.tbl['-CO'])
    } else if(!'hm.color' %in% existing.vl){
        updated.vl$hm.color <- "default"
    }

    #### Color distribution. ####
    if('-CD' %in% names(args.tbl)) {
        stopifnot(as.numeric(args.tbl['-CD']) > 0)
        updated.vl$color.distr <- as.numeric(args.tbl['-CD'])
    } else if(!'color.distr' %in% existing.vl) {
        updated.vl$color.distr <- .6
    }

    #### Low count cutoff for rank-based normalization ####
    if('-LOW' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-LOW']) >= 0)
        updated.vl$low.count <- as.integer(args.tbl['-LOW'])
    } else if(!'low.count' %in% existing.vl) {
        updated.vl$low.count <- 10
    } else {  # ensure low.count is not empty.
        updated.vl$low.count <- low.count
    }
    if(!is.null(low.count)) {
        updated.vl$low.count.ratio <- updated.vl$low.count / low.count
    }

    #### Misc. options for heatmap. ####
    updated.vl$go.paras <- go.paras
    if('-KNC' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-KNC']) > 0)
        updated.vl$go.paras$knc <- as.integer(args.tbl['-KNC'])
    } else if(!'knc' %in% names(go.paras)) {
        updated.vl$go.paras$knc <- 5
    }
    if('-MIT' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-MIT']) > 0)
        updated.vl$go.paras$max.iter <- as.integer(args.tbl['-MIT'])
    } else if(!'max.iter' %in% names(go.paras)) {
        updated.vl$go.paras$max.iter <- 20
    }
    if('-NRS' %in% names(args.tbl)) {
        stopifnot(as.integer(args.tbl['-NRS']) > 0)
        updated.vl$go.paras$nrs <- as.integer(args.tbl['-NRS'])
    } else if(!'nrs' %in% names(go.paras)) {
        updated.vl$go.paras$nrs <- 30
    }


    updated.vl
}

CheckHMColorConfig <- function(hm.color, bam.pair) {
    if(hm.color != 'default') {
        v.colors <- unlist(strsplit(hm.color, ":"))
        if(bam.pair && length(v.colors) != 2 && length(v.colors) != 3 || 
           !bam.pair && length(v.colors) != 1 && length(v.colors) != 2) {
            stop("Heatmap color specifications must correspond to bam-pair!\n")
        }
    }
}

EchoPlotArgs <- function() {
    cat("## Plotting-related parameters:\n")
    cat("### Misc. parameters:\n")
    cat("    -FS Font size(default=12)\n")
    cat("### Avg. profiles parameters:\n")
    cat("    -WD Image width(default=8)\n")
    cat("    -HG Image height(default=7)\n")
    cat("    -SE  Shall standard errors be plotted?(0 or 1)\n")
    cat("    -MW  Moving window width to smooth avg. profiles, must be integer\n")
    cat("           1=no(default); 3=slightly; 5=somewhat; 9=quite; 13=super.\n")
    cat("    -H   Opacity of shaded area, suggested value:[0, 0.5]\n")
    cat("           default=0, i.e. no shading, just lines\n")
    cat("    -YAS Y-axis scale: auto(default) or min_val,max_val(custom scale)\n")
    cat("    -LEG Draw legend? 1(default) or 0\n")
    cat("    -BOX Draw box around plot? 1(default) or 0\n")
    cat("    -VLN Draw vertical lines? 1(default) or 0\n")
    cat("    -XYL Draw X- and Y-axis labels? 1(default) or 0\n")
    cat("    -LWD Line width(default=3)\n")
    cat("### Heatmap parameters:\n")
    cat("    -GO  Gene order algorithm used in heatmaps: total(default), hc, max,\n")
    cat("           prod, diff, km and none(according to gene list supplied)\n")
    cat("    -LOW Low count cutoff(default=10) in rank-based normalization\n")
    cat("    -KNC K-means or HC number of clusters(default=5)\n")
    cat("    -MIT Maximum number of iterations(default=20) for K-means\n")
    cat("    -NRS Number of random starts(default=30) in K-means\n")
    cat("    -RR  Reduce ratio(default=30). The parameter controls the heatmap height\n")
    cat("           The smaller the value, the taller the heatmap\n")
    cat("    -SC  Color scale used to map values to colors in a heatmap\n")
    cat("           local(default): base on each individual heatmap\n")
    cat("           region: base on all heatmaps belong to the same region\n")
    cat("           global: base on all heatmaps together\n")
    cat("           min_val,max_val: custom scale using a pair of numerics\n")
    cat("    -FC  Flooding fraction:[0, 1), default=0.02\n")
    cat("    -CO  Color for heatmap. For bam-pair, use color-tri(neg_color:[neu_color]:pos_color)\n")
    cat("           Hint: must use R colors, such as darkgreen, yellow and blue2\n")
    cat("                 The neutral color is optional(default=black)\n")
    cat("    -CD  Color distribution for heatmap(default=0.6). Must be a positive number\n")
    cat("           Hint: lower values give more widely spaced colors at the negative end\n")
    cat("                 In other words, they shift the neutral color to positive values\n")
    cat("                 If set to 1, the neutral color represents 0(i.e. no bias)\n")

}

EchoCoverageArgs <- function() {
    cat("## Coverage-generation parameters:\n")
    cat("  -F   Further information provided to select database table or plottype:\n")
    cat("         This is a string of description separated by comma.\n")
    cat("         E.g. protein_coding,K562,rnaseq(order of descriptors does not matter)\n")
    cat("              means coding genes in K562 cell line drawn in rnaseq mode.\n")
    cat("  -D   Gene database: ensembl(default), refseq\n")
    cat("  -L   Flanking region size(will override flanking factor)\n")
    cat("  -N   Flanking region factor\n")
    cat("  -RB  The fraction of extreme values to be trimmed on both ends\n")
    cat("         default=0, 0.05 means 5% of extreme values will be trimmed\n")
    cat("  -S   Randomly sample the regions for plot, must be:(0, 1]\n")
    cat("  -P   #CPUs to use. Set 0(default) for auto detection\n")
    cat("  -AL  Algorithm used to normalize coverage vectors: spline(default), bin\n")
    cat("  -CS  Chunk size for loading genes in batch(default=100)\n")
    cat("  -MQ  Mapping quality cutoff to filter reads(default=20)\n")
    cat("  -FL  Fragment length used to calculate physical coverage(default=150)\n")
    cat("  -SS  Strand-specific coverage calculation: both(default), same, opposite\n")
    cat("  -IN  Shall interval be larger than flanking in plot?(0 or 1, default=automatic)\n")
    cat("  -FI  Forbid image output if set to 1(default=0)\n")
}
############### End arguments configuration #####################
#################################################################
