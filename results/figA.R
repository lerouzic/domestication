#!/usr/bin/env Rscript

# First figure: dynamics of the effective pop size and average fitness +/- sd for a random population

source("../src/makeparam_functions.R")

out.dir <- "../cache/simDefault"
target.rep <- NA # If NA, pick a random replicate

get.N.dyn <- function(repdir) {
	ll <- list.files(path=repdir, pattern="^param.*\\.txt$", full.names=TRUE)
	myN <- sapply(ll, function(ff) { 
		if (length(grep("INIT_PSIZE", readLines(ff))) == 0) 
			return(NA)
		else {
			pp <- read.param(ff)
			return(pp[["INIT_PSIZE"]])
		}
		})
	myN <- na.omit(myN)[cumsum(!is.na(myN))] # Replace NAs by the last non-NA
	names(myN) <- regmatches(ll, regexpr(ll, pattern="gen(\\d+)"))
	myN
}



if (is.na(target.rep)) {
	rep.dirs <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	if (length(rep.dirs) == 0) stop("Unable to find the results in ", out.dir, ".")
	target.rep <- sample(rep.dirs, 1)
}

out.files <- list.files(pattern="out.*", path=target.rep, full.names=TRUE)

stopifnot(length(out.files)>0)

out.mean <- replicate.apply(out.files, mean)
out.var  <- replicate.apply(out.files, var)

mfit <- out[,"MFit"]
vfit <- out[,"VFit"]
gen <- as.numeric(out[,"Gen"])

N <- get.N.dyn(target.rep)[gen+1]

plot(gen, N, type="l", xlim=range(gen), ylim=c(0, max(N)))
lines(gen, N/(1+4*vfit/mfit/mfit), col="blue")


plot(NULL, xlab="Generations", ylab="Fitness", xlim=range(gen), ylim=c(0,1))
lines(gen, mfit+sqrt(vfit), lty=1, col="blue")
lines(gen, mfit-sqrt(vfit), lty=1, col="blue")
lines(gen, mfit, lwd=2)

