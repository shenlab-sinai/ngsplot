#!/usr/bin/env Rscript
# Install the newly obtained DB tables into two ngs.plot system files:
# default.tbl and dbfile.tbl.

args <- commandArgs(T)
# args <- c('/tmp/tmpLXAGXl', '/tmp/tmpUJduMe')
default.tbl <- args[1]
db.list <- args[2]

# Genome-region default values.
anno.tbl <- read.table(default.tbl, header=TRUE, sep="\t",
                       blank.lines.skip=TRUE, stringsAsFactors=FALSE)
row.names(anno.tbl) <- paste(anno.tbl$Genome, anno.tbl$Region, sep=".")

# Database table list.
anno.db.tbl <- read.table(db.list, 
                          col.names=c("db.file", "Genome", "DB", "Region", 
                                      "FI.1", "FI.2", "FI.3"),
                          sep="\t", blank.lines.skip=TRUE, 
                          stringsAsFactors=FALSE)
tss.tbl <- anno.db.tbl[anno.db.tbl$Region=="genebody", ]
tss.tbl$Region <- "tss"
tes.tbl <- anno.db.tbl[anno.db.tbl$Region=="genebody", ]
tes.tbl$Region <- "tes"
anno.db.tbl <- rbind(anno.db.tbl, tss.tbl, tes.tbl)
anno.db.tbl$URL <- NA  # reserve for future use.

# Extract further info columns.
default.fi <- anno.tbl[, c("DefaultFI1", "DefaultFI2", "DefaultFI3")]
db.fi <- anno.db.tbl[, c("FI.1", "FI.2", "FI.3")]
db.fi$tag <- paste(anno.db.tbl$Genome, anno.db.tbl$Region, sep=".")

# DB scores.
getScore <- function(v.fi) {
# Calculate DB score based on intersection with default values.
# Args:
#   v.fi: character vector of further infos.
# Returns: integer of DB score.

    default.vals <- default.fi[v.fi["tag"], ]
    fi.vals <- v.fi[c("FI.1", "FI.2", "FI.3")]
    3 - length(intersect(default.vals, fi.vals))
}
db.score <- apply(db.fi, 1, getScore)
anno.db.tbl$dbScore <- db.score

# Save to text output.
progpath <- Sys.getenv('NGSPLOT')
if(progpath == "") {
	stop("Set environment variable NGSPLOT before run the program. See README for details.\n")
}
def.f <- file.path(progpath, 'database', 'default.tbl')
db.f <- file.path(progpath, 'database', 'dbfile.tbl')
org.anno.tbl <- read.delim(def.f)
org.anno.db.tbl <- read.delim(db.f)
anno.tbl <- rbind(org.anno.tbl, anno.tbl)
anno.tbl <- unique(anno.tbl)
anno.db.tbl <- rbind(org.anno.db.tbl, anno.db.tbl)
anno.db.tbl <- unique(anno.db.tbl)
write.table(anno.tbl, file=def.f, row.names=F, sep="\t", quote=F)
write.table(anno.db.tbl, file=db.f, row.names=F, sep="\t", quote=F)

# getScore <- function(i, db.fi, default.fi){
#     score <- 3 - length(intersect(as.matrix(default.fi[db.fi[i, "tag"], ]), as.matrix(db.fi[i, ])))
#     return(score)
# }

# for (i in (1:dim(db.fi)[1])){
#     anno.db.tbl[i, "dbScore"] <- getScore(i, db.fi, default.fi)
# }

# anno.version="0.01"
# save(anno.tbl, anno.db.tbl, anno.version, file="database.RData")
