# Various functions to analyze the network


source("../src/analysis_tools.R")

library(Rcpp)
library(inline, quietly=TRUE)

cppFunction('
	List internal_loop_cpp(const NumericMatrix &W, const NumericVector &S0, double a, double env, unsigned int steps, unsigned int measure) {
		double lambda = (1-a)/a;
		double mu     = 1/(a*(1-a));
		NumericMatrix sto (S0.size(), steps+1);
		NumericVector sumx (S0.size());
		NumericVector sumx2 (S0.size()); 
		sto(0,0) = env;
		for (unsigned int i = 1; i < S0.size(); i++)
			sto(i,0) = S0(i);
		for (unsigned int t = 1; t <= steps; t++) {
			for (unsigned int i = 0; i < S0.size(); i++) {
				double tmp = 0.;
				if (i == 0) { // The signal gene
					tmp = S0[i];
				} else {
					for (unsigned int j = 0; j < S0.size(); j++) {
						tmp += sto(j,t-1) * W(i,j);
					}
					tmp =  1. / (1. + lambda * exp(-mu*tmp));
				}
				sto(i,t) = tmp;
				if (t > steps-measure) {
					sumx(i) += tmp;
					sumx2(i) += tmp*tmp;
				}
			}
		}
		for (unsigned int i = 0; i < S0.size(); i++) {
			sumx(i) /= static_cast<double>(measure); // sumx(i) now contains the mean
			sumx2(i) /= static_cast<double>(measure);
			sumx2(i) += -sumx(i)*sumx(i); // sumx2(i) now contains the variance
		}
		return List::create(Named("full")=sto, Named("mean")=sumx, Named("var")=sumx2);
	}')


model.M2 <- function(W, S0=rep(a, nrow(W)), a=0.2, env=0.5, steps=20, measure=4, full=FALSE) {
    ans <- internal_loop_cpp(W, S0, a, env, steps, measure)
	if (!full) ans$full <- NULL
	return(ans)
}

cleanW <- function(W, epsilon=NULL, env=0.5, ...) {
	if (is.null(epsilon)) epsilon <- sqrt((nrow(W)-1)*0.01^2)
	cleanW <- matrix(0, nrow=nrow(W), ncol=ncol(W))
	ref.expr <- model.M2(W, env=env, ...)$mean
	for (i in 1:nrow(W))
		for (j in 1:ncol(W)) {
			myW <- W
			myW[i,j] <- 0
			new.expr <- model.M2(myW, env=env, ...)$mean
			dd <- sqrt(sum((ref.expr[-1]-new.expr[-1])^2))
			if (dd > epsilon)
				cleanW[i,j] <- W[i,j]
		}
	cleanW
}

number.connections <- function(W, epsilon=NULL, env=0.5, ...) {
	cW <- cleanW(W=W, epsilon=epsilon, env=env, ...)
	sum(cW != 0)
}

number.connections.dyn <- function(out.table, epsilon=NULL, env=0.5) { # if env == NULL, MPhen1 is used instead
	W.table <- out.table[,grepl(colnames(out.table), pattern="MeanAll")]
	env <- if(is.null(env)) out.table[,"MPhen1"] else rep(env, nrow(out.table))
	
	net.size <- sqrt(ncol(W.table))
	nb.conn <- sapply(1:nrow(W.table), function(i) { 
			W <- matrix(unlist(W.table[i,]), ncol=net.size, byrow=TRUE)
			stopifnot(nrow(W) == ncol(W))
			number.connections(W, epsilon=epsilon, env=env[i])
		})
	names(nb.conn) <- as.character(out.table[,"Gen"])
	nb.conn
}

# Average out all network connections from a directory 
mean.connect <- function(out.dir, env=0.5, epsilon=NULL, max.reps=Inf, mc.cores=detectCores()-1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)

	tt <- mclapply(out.files[1:(min(max.reps, length(out.files)))], read.table, header=TRUE, mc.cores=min(mc.cores, 4)) 
	nn <- mclapply(tt, number.connections.dyn, env=env, epsilon=epsilon, mc.cores=mc.cores)
	ans <- colMeans(do.call(rbind, nn))
	rm(tt)
	gc()
	return(ans)
}

