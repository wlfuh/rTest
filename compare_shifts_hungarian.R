compass_score_hungry <- function (P, Q) {
	require(pracma)
	require(clue)   
	# computes the COMPASS score between points in 2D space 
	stopifnot(is.numeric(P), is.numeric(Q))
	if (is.vector(P)) 
			P <- matrix(P, ncol = 1)
	if (is.vector(Q)) 
			Q <- matrix(Q, ncol = 1)
	if (ncol(P) != ncol(Q)) 
			stop("'P' and 'Q' must have the same number of columns.")
	D <- distmat(P, Q)
	assignments <- solve_LSAP(D)				
	return(mean(D[cbind(seq_along(assignments), assignments)]))
}

load("correlation_TOCSY.RData")
haus <- NULL
P <- as.matrix(cs2d_exp[,c("cs1","cs2")])
frames <- unique(cs2d_pred$model)
for (i in seq_along(frames)){
	frame <- frames[i]
	cs2d_pred_tmp <- subset(cs2d_pred,model==frame)
	Q <- as.matrix(cs2d_pred_tmp[,c("cs1","cs2")])
	haus <- c(haus,compass_score_hungry(P,Q))
	cat(sprintf("frame: %s COMPASS: %s\n",frame,haus[i]))
}
save.image("compass_score_hungry_comparison.RData")
