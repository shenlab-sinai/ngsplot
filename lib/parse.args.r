# Parse command line arguments and extract them into an associate array.
# Check if the required arguments are all satisfied.
parse.args <- function(args, manditories){
	if(length(args) %% 2 == 1 || length(args) == 0){
		warning('Unpaired argument and value.')
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
		warning(paste('Missing argument: ', paste(manditories[!man.bool], collapse=','), '.', sep=''))
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
