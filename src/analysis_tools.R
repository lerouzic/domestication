# Various functions to analyse the simulation results

# Applies the function FUN (typically, mean or var) to each element of the matrix/data.frame stored in a vector of files
replicate.apply <- function(files, FUN=mean, ...) {
	library(abind)
	tt <- lapply(files, function(ff) read.table, header=TRUE)
	# better checking if the data base is consistent
	stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
	arr <- do.call(abind, tt, along=3)
	apply(arr, 1:2, FUN)
}

# Reconstructs the dynamics of the population size N from a replicate directory 
# (takes a bit of time due to the number of files)
get.Ndyn <- function(repdir) {
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

# Detects the begining and the end of the bottleneck(s?) from the dynamics of N
bottleneck.detect <- function(Ndyn) {
	# Returns the generation number if Ndyn is named with generations, otherwise the index
	if (is.null(names(Ndyn)) names(Ndyn) <- as.character(seq_along(Ndyn))
	stopifnot(length(Ndyn) < 3)
	if (length(unique(Ndyn)) > 10) warning("Probably too many changes in N to detect a bottleneck")
	bottle.begin <- which(diff(Ndyn) < 0)
	bottle.end <- which(diff(Ndyn) > 0)
	list(begin=as.numeric(names(Ndyn)[bottle.begin]), end=as.numeric(names(Ndyn)[bottle.end]))
}

# Tries to estimate the selection regime for each gene. The output code is s for stabilizing, p for plastic, n for neutral (non-selected). In addition, c is 
# used for constant (i.e. stable genes which optimum does not change during domestication). 
# The output is two characters, one before and one after domestication (even for constant genes, which are cc for consistency). 
selectionregime.detect <- function(out.table) {
	# Assumes only one change in the selection regime
	opts <- out.table[,grepl(rownames(out.table), pattern="theta")]
	stopifnot (length(opts) > 0)
	
	stability <- lapply(opts function(x) { ans <- ifelse(diff(x) == 0, "s", "p"); ans[x[-1] == 0] <- "n"; paste(rle(ans)$values, collapse="") } 
	category <- sapply(stability, function(x) if (x == "s") "cc" else if (x == "p") "pp" else if (x == "n") "n" else x)
	return(category)
}

# Detects the generation number(s?) at which the selection regime changes
selectionchange.detect <- function(out.table) {
	# The task is difficult, has to rely on some assumptions 
	# The algo looks for stable optima that changes from time to time
	# out.file can be the mean over replicates or a single file. Probably works better with a single file
	
	selcat <- selectionregime.detect(out.table)
	changes.domestic <- lapply(opts[selcat == ss], function(x) which(diff(x)!=0))
	if (!do.call(all.equal, changes.domestic)) stop("The selection change history looks different for domestic genes.")
	return(out.table$Gen[changes.domestic[[1]]])
}

# Adds the bottleneck mark in a figure
bottleneck.plot <- function(Ndyn, y=0, code=3, angle=90, ...) {
	bttl <- bottleneck.detect(Ndyn)
	arrows(x0=bbtl$begin, x1=bbtl$end, y0=y, code=code, angle=angle, ...)
}

# Adds the selection mark in a figure
selectionchange.plot <- function(out.table, y=0, pch=25, col="black", bg="red", ...) {
	sel <- selectionchange.detect(out.table)
	points(x=sel, y=y, pch=pch, col=col, bg=bg, ...)
}
