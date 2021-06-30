# Variance in relative fitness


relvfit <- function(out.table) {
	out.table[,"VFit"]/(out.table[,"MFit"]^2)
}

relvfit.dir <- function(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern=NULL) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	tt <- results.table(out.files, mc.cores, max.reps, colnames.pattern=colnames.pattern)
	ans <- lapply(tt, relvfit)
	rm(tt)
	gc()
	return(ans)
}

relvfit.mean <- function(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern=NULL) {
	data <- relvfit.dir(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern=NULL)
	colMeans(do.call(rbind, data))
}

relvfit.mean.cache <- function(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern=NULL)
	cache.fun(relvfit.mean, out.dir=out.dir, max.reps=max.reps, mc.cores=mc.cores, colnames.pattern=colnames.pattern, cache.subdir="Rcache-relvfit", file.prefix=basename(out.dir))

relvfit.quantile <- function(out.dir, quant=0.5, max.reps=Inf, mc.cores=1, colnames.pattern=NULL) {
	data <- relvfit.dir(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern=NULL)
	apply(do.call(rbind, data), 2, quantile, prob=quant)
}

relvfit.quantile.cache <- function(out.dir, quant=0.5, max.reps=Inf, mc.cores=1, colnames.pattern=NULL)
	cache.fun(relvfit.quantile, out.dir=out.dir, quant=quant, max.reps=max.reps, mc.cores=mc.cores, colnames.pattern=colnames.pattern, cache.subdir="Rcache-relvfit", file.prefix=basename(out.dir) )
