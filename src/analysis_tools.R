# Various functions to analyse the simulation results

suppressMessages(library(abind))

library(parallel)

source("../src/makeparam_functions.R") # dubious path management...
source("../src/cache.R")

out.files <- function(out.dir) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- paste0(out.reps, "/out-", sub(basename(out.reps), pattern="-", replacement=""), ".txt")
	out.files[file.exists(out.files)]
}

# Moving average (generation numbers as names)
mov.avg <- function(x, gen=1:length(x), size=1, min.gen=0) {
	x.ma <- diff(cumsum(x), lag=size)/size
	gen.ma <- diff(cumsum(gen), lag=size)/size
	names(x.ma) <- as.character(round(gen.ma))
	x.ma <- x.ma[as.numeric(x.ma) >= min.gen]
	x.ma <- x.ma[!is.na(x.ma)]
	x.ma
}

# Applies the function FUN (typically, mean or var) to each element of the list of matrix/data.frames 
replicate.apply <- function(tt, FUN=mean, ...) {
	# better checking if the data base is consistent
	stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
	arr <- do.call(abind, c(tt, list(along=3)))
	ans <- apply(arr, 1:2, FUN, ...)
	ans
}

# This is faster than replicate.apply(tt, mean)
replicate.mean <- function(tt) {
	stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
	arr <- do.call(abind, c(tt, list(along=3)))
	rowMeans(arr, dims=2)
}

# This is faster than replicate.apply(tt, var)
replicate.var <- function(tt) {
	# This rowVars function comes from https://stat.ethz.ch/pipermail/r-help/2006-April/103001.html
	# Author: David Brahm
	.rowVars <- function(x, na.rm=FALSE, dims=1, unbiased=TRUE, SumSquares=FALSE, twopass=FALSE) {
		if (SumSquares) return(rowSums(x^2, na.rm, dims))
		N <- rowSums(!is.na(x), FALSE, dims)
		Nm1 <- if (unbiased) N-1 else N
		if (twopass) {x <- if (dims==0) x - mean(x, na.rm=na.rm) else sweep(x, 1:dims, rowMeans(x,na.rm,dims))}
		(rowSums(x^2, na.rm, dims) - rowSums(x, na.rm, dims)^2/N) / Nm1
	}
	stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
	arr <- do.call(abind, c(tt, list(along=3)))
	.rowVars(arr, dims=2)
}

replicate.quantile <- function(tt, quant) {
	replicate.apply(tt, quantile, prob=quant)
}

# Reconstructs the dynamics of the population size N from a replicate directory 
# (takes a bit of time due to the number of files)
get.Ndyn <- function(repdir) {
	ll <- list.files(path=repdir, pattern="^param.*\\.txt$", full.names=TRUE)
	cc <- numeric(0)
	if (length(ll) == 0) {
		cc <- list.files(path=repdir, pattern="^param.*\\.tar", full.names=TRUE)
		stopifnot(length(cc) == 1)
		untar.param(cc, remove.tar=FALSE)
		ll <- list.files(path=repdir, pattern="^param.*\\.txt$", full.names=TRUE)
	}
	myN <- sapply(ll, function(ff) { 
		if (length(grep("INIT_PSIZE", readLines(ff))) == 0) 
			return(NA)
		else {
			pp <- read.param(ff)
			return(pp[["INIT_PSIZE"]])
		}
		})
	if (length(cc) == 1) # there was a tar that was uncompressed
		unlink(ll)
	myN <- na.omit(myN)[cumsum(!is.na(myN))] # Replace NAs by the last non-NA
	names(myN) <- sapply(regmatches(ll, regexec(ll, pattern="gen(\\d+)")), "[", 2)
	names(myN) <- as.character(as.numeric(names(myN)))
	myN
}

get.Ndyn.cache <- function(repdir) {
	cache.fun(get.Ndyn, repdir=repdir, cache.subdir="Rcache-Ndyn", file.prefix=basename(repdir))
}

# Detects the begining and the end of the bottleneck(s?) from the dynamics of N
bottleneck.detect <- function(Ndyn) {
	# Returns the generation number if Ndyn is named with generations, otherwise the index
	if (is.null(names(Ndyn))) names(Ndyn) <- as.character(seq_along(Ndyn))
	stopifnot(length(Ndyn) > 3)
	if (length(unique(Ndyn)) > 10) warning("Probably too many changes in N to detect a bottleneck")
	bottle.begin <- which(diff(Ndyn) < 0) + 1
	bottle.end <- which(diff(Ndyn) > 0) + 1
	list(begin=as.numeric(names(Ndyn)[bottle.begin]), end=as.numeric(names(Ndyn)[bottle.end]))
}

# Tries to estimate the selection regime for each gene. The output code is s for stabilizing, p for plastic, n for neutral (non-selected). In addition, c is 
# used for constant (i.e. stable genes which optimum does not change during domestication). 
# The output is two characters, one before and one after domestication (even for constant genes, which are cc for consistency). 
selectionregime.detect <- function(out.table) {
	# Assumes only one change max in the selection regime
	opts <- out.table[,grepl(colnames(out.table), pattern="FitOpt")]
	stopifnot (length(opts) > 0)
	stability <- lapply(as.data.frame(opts), function(x) { ans <- ifelse(diff(x) == 0, "s", "p"); ans[x[-1] == 0] <- "n"; paste(rle(ans)$values, collapse="") } )
	if (all(nchar(stability)==1)) { # No change in selection regime
		category <- unlist(stability)
	} else {
		category <- sapply(stability, function(x) if (x == "s") "cc" else if (x == "p") "pp" else if (x == "n") "nn" else if (x=="nps") "ns" else if (x == "sps") "ss" else x)
	}
	return(category)
}

# Detects the generation number(s?) at which the selection regime changes
selectionchange.detect <- function(out.table) {
	# The task is difficult, has to rely on some assumptions 
	# The algo looks for stable optima that changes from time to time
	# out.file can be the mean over replicates or a single file. Probably works better with a single file
	selcat <- selectionregime.detect(out.table)
	changes.domestic <- as.matrix(apply(out.table[,names(selcat)[selcat == "ss"],drop=F], 2, function(x) which(diff(x)!=0)))
	if (!all(apply(changes.domestic, 2, function(x) length(unique(x))==1))) stop("The selection change history looks different for domestic genes.")
	return(out.table[1+changes.domestic[1,],"Gen"])
}

# Adds the bottleneck mark in a figure
bottleneck.plot <- function(Ndyn, y=0, code=3, angle=90, length=0.05, ...) {
	bttl <- bottleneck.detect(Ndyn)
	arrows(x0=bttl$begin, x1=bttl$end, y0=y, code=code, angle=angle, length=length, ...)
}

# Adds the selection mark in a figure
selectionchange.plot <- function(out.table, y=0, pch=25, col="black", bg="red", ...) {
	sel <- selectionchange.detect(out.table)
	points(x=sel, y=y, pch=pch, col=col, bg=bg, ...)
}

# Returns a list of data.frames from all files (without truncated files)
results.table <- function(out.files, mc.cores=1, max.reps=length(out.files), verbose=TRUE, colnames.pattern=NULL) {
	tt <- mclapply(out.files[1:(min(max.reps, length(out.files)))], function(ff) { 
		ans <- try(read.table(ff, header=TRUE))
		if (class(ans) == "try-error") return(numeric(0))
		if ("Gen" %in% colnames(ans))
			rownames(ans) <- as.character(ans$Gen)
		if (!is.null(colnames.pattern))
			ans <- ans[,grepl(pattern=colnames.pattern, colnames(ans))]
		return(ans)
		}, mc.cores=mc.cores)
	ltt <- sapply(tt, nrow)
	full.length <- max(ltt)
	if (verbose && any(ltt != full.length))
		cat("\nOnly ", sum(ltt == full.length), " out of ", length(out.files), " were included in the analysis\n")
	tt <- tt[ltt == full.length]
	tt
}

# Average out all data tables from a directory
mean.sim <- function(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern=NULL) {
	out.files <- out.files(out.dir)
	tt <- results.table(out.files, mc.cores, max.reps, colnames.pattern=colnames.pattern)
	ans <- replicate.mean(tt)
	rm(tt)
	gc()
	return(ans)
}

mean.sim.cache <- function(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern=colnames.pattern) {
	cache.fun(mean.sim, out.dir=out.dir, max.reps=max.reps, mc.cores=mc.cores, colnames.pattern=colnames.pattern, cache.subdir="Rcache-means", file.prefix=basename(out.dir))
}

quantile.sim <- function(out.dir, quant=0.5, max.reps=Inf, mc.cores=1, colnames.pattern=NULL) {
	out.files <- out.files(out.dir)
	tt <- results.table(out.files, mc.cores, max.reps, colnames.pattern=colnames.pattern)
	ans <- replicate.quantile(tt, quant)
	rm(tt)
	gc()
	return(ans)
}

quantile.sim.cache <- function(out.dir, quant=0.5, max.reps=Inf, mc.cores=1, colnames.pattern=NULL) {
	cache.fun(quantile.sim, out.dir=out.dir, quant=quant, max.reps=max.reps, mc.cores=mc.cores, colnames.pattern=colnames.pattern, cache.subdir="Rcache-quantiles", file.prefix=basename(out.dir))
}


# variance of all data tables from a directory
var.sim <- function(out.dir, max.reps=Inf, mc.cores=detectCores()-1) {
	out.files <- out.files(out.dir)
	tt <- results.table(out.files, mc.cores, max.reps)
	ans <- replicate.var(tt)
	rm(tt)
	gc()
	return(ans)
}

