#!/usr/bin/env Rscript

res.time <- system.time({  # Load required libraries.
    if(!suppressMessages(require(ShortRead, warn.conflicts=F))) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(ShortRead)
        if(!suppressMessages(require(ShortRead, warn.conflicts=F))) {
            stop('Loading package ShortRead failed!')
        }
    }
    if(!suppressMessages(require(rtracklayer, warn.conflicts=F))) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(rtracklayer)
        if(!suppressMessages(require(rtracklayer, warn.conflicts=F))) {
            stop('Loading package rtracklayer failed!')
        }
    }
    if(!suppressMessages(require(BSgenome, warn.conflicts=F))) {
        if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
        BiocManager::install(BSgenome)
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
})

print(res.time)
