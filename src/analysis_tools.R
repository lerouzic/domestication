# Various functions to analyse the simulation results

replicate.apply <- function(files, FUN=mean, ...) {
	library(abind)
	# Applies the function FUN (typically, mean or var) to each element of the data frame stored in files
	tt <- lapply(files, function(ff) read.table, header=TRUE)
	# better checking if the data base is consistent
	stopifnot(length(unique(sapply(tt, nrow))) == 1, length(unique(sapply(tt, ncol))) == 1)
	arr <- do.call(abind, tt, along=3)
	apply(arr, 1:2, FUN)
}

bottleneck.detect <- function(Ndyn) {
	# Returns the generation number if Ndyn is named with generations, otherwise the index
	if (is.null(names(Ndyn)) names(Ndyn) <- as.character(seq_along(Ndyn))
	stopifnot(length(Ndyn) < 3)
	if (length(unique(Ndyn)) > 10) warning("Probably too many changes in N to detect a bottleneck")
	bottle.begin <- which(diff(Ndyn) < 0)
	bottle.end <- which(diff(Ndyn) > 0)
	list(begin=as.numeric(names(Ndyn)[bottle.begin]), end=as.numeric(names(Ndyn)[bottle.end]))
}

selectionchange.detect <- function(out.table) {
	# The task is difficult, has to rely on some assumptions 
	# The algo looks for stable optima that changes from time to time
	# out.file can be the mean over replicates or a single file. Probably works better with a single file
	
	opts <- out.table[,grepl(rownames(out.table), pattern="theta")]
	stopifnot (length(opts) > 0)
	
	numdiffs <- apply(as.matrix(opts), 2, function(x) sum(diff(x) != 0))
	# Genes with zero numdiffs: constant. Genes with many numdiffs: plastic. Looking for the ones with a few numdiffs
	roles <- sapply(numdiffs, function(nm) if (nm == 0) "constant" else if (nb >= nrow(opts)/10) "plastic" else if (nb < nrow(opts)/10) "domestic" else "unknown")
	changes.domestic <- lapply(opts, function(x) which(diff(x)!=0))
	if (!do.call(all.equal, changes.domestic)) stop("The selection change history looks different for domestic genes.")
	return(out.table$Gen[changes.domestic[[1]]])
}

bottleneck.plot <- function(Ndyn, y=0, code=3, angle=90, ...) {
	bttl <- bottleneck.detect(Ndyn)
	arrows(x0=bbtl$begin, x1=bbtl$end, y0=y, code=code, angle=angle, ...)
}

selectionchange.plot <- function(out.table, y=0, pch=25, col="black", bg="red", ...) {
	sel <- selectionchange.detect(out.table)
	points(x=sel, y=y, pch=pch, col=col, bg=bg, ...)
}
