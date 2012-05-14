# Separate a file name into directory and base file name.
sep.filename <- function(f){
	p <- unlist(strsplit(f, '/'))
	if(length(p) == 1){
		dir <- '.'
		n <- p
	}else{
		dir <- paste(p[1:(length(p)-1)], collapse='/')
		n <- p[length(p)]
	}
	list(path=dir, name=paste(n,"$",sep=""))
}
