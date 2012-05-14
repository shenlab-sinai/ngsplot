# Smooth plotting data in a vector.
smoothvec <- function(v, radius, method=c('mean', 'median')){
	stopifnot(is.vector(v))
	stopifnot(length(v) > 0)
	stopifnot(radius > 0 && radius < 1)
	halfwin <- ceiling(length(v) * radius)
	s <- vector()
	for(i in 1:length(v)){
		winpos <- (i - halfwin) : (i + halfwin)
		winpos <- winpos[winpos > 0 & winpos <= length(v)]
		if(method == 'mean'){
			s[i] <- mean(v[winpos])
		}else if(method == 'median'){
			s[i] <- median(v[winpos])
		}
	}
	s
}

# Smooth plotting data which can be a matrix or a vector.
smoothplot <- function(m, radius, method=c('mean', 'median')){
	stopifnot(is.matrix(m) || is.vector(m))
	if(is.matrix(m)){
		for(i in 1:ncol(m)){
			m[, i] <- smoothvec(m[, i], radius, method)
		}
	}else{
		m <- smoothvec(m, radius, method)
	}
	m
}
