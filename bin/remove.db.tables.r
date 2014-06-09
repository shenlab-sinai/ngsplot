#!/usr/bin/env Rscript
# Remove the genome entries from two ngs.plot system files:
# default.tbl and dbfile.tbl.

args <- commandArgs(T)
if(length(args)>1){
	gname.to.rm <- args[1]
	feature.to.rm <- args[2]
}else{
	gname.to.rm <- args[1]
}

# Save to text output.
progpath <- Sys.getenv('NGSPLOT')
if(progpath == "") {
	stop("Set environment variable NGSPLOT before run the program. See README for details.\n")
}
default.file <- file.path(progpath, 'database', 'default.tbl')
dbfile.file <- file.path(progpath, 'database', 'dbfile.tbl')
default.tbl <- read.delim(default.file)
dbfile.tbl <- read.delim(dbfile.file)

if(exists("feature.to.rm")){
	default.tbl <- subset(default.tbl, !(Genome==gname.to.rm & Region==feature.to.rm))
	dbfile.tbl <- subset(dbfile.tbl, !(Genome==gname.to.rm & Region==feature.to.rm))
}else{
	default.tbl <- default.tbl[default.tbl$Genome != gname.to.rm, ]
	dbfile.tbl <- dbfile.tbl[dbfile.tbl$Genome != gname.to.rm, ]	
}


write.table(default.tbl, file=default.file, row.names=F, sep="\t", quote=F)
write.table(dbfile.tbl, file=dbfile.file, row.names=F, sep="\t", quote=F)


# getScore <- function(i, db.fi, default.fi){
#     score <- 3 - length(intersect(as.matrix(default.fi[db.fi[i, "tag"], ]), as.matrix(db.fi[i, ])))
#     return(score)
# }

# for (i in (1:dim(db.fi)[1])){
#     anno.db.tbl[i, "dbScore"] <- getScore(i, db.fi, default.fi)
# }

# anno.version="0.01"
# save(anno.tbl, anno.db.tbl, anno.version, file="database.RData")
