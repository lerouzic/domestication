# Network fearures: applies a function to the W matrix

WFUN.dyn <- function(out.table, WFUN="mean", deltaG=1) {
	gen <- out.table[,"Gen"]
	seqgen <- seq(1, length(gen), length.out=min(length(gen), 1+gen[length(gen)] %/% deltaG))
	listW <- Wlist.table(out.table[seqgen,])
	ww <- unlist(lapply(listW, eval(parse(text=WFUN))))
	ww
}

WFUN.dir <- function(out.dir, WFUN="mean", deltaG=1, mc.cores=1) {
	out.files <- out.files(out.dir)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		wm <- WFUN.dyn(tt, WFUN=WFUN, deltaG=deltaG)
		rm(tt)
		wm
	}, mc.cores=mc.cores)
	ansl <- sapply(ans, length)
	ans <- ans[ansl==max(ansl)]
	ans
}

WFUN.mean <- function(out.dir, WFUN="mean", deltaG=1, mc.cores=1) {
	data <- WFUN.dir(out.dir, WFUN=WFUN, deltaG=deltaG, mc.cores=mc.cores)
	colMeans(do.call(rbind, data))
}

WFUN.mean.cache <- function(out.dir, WFUN="mean", deltaG=1, mc.cores=1, cache.subdir=paste0("Rcache-W", WFUN)) {
	cache.fun(WFUN.mean, out.dir=out.dir, WFUN=WFUN, deltaG=deltaG, mc.cores=mc.cores, cache.subdir=cache.subdir, file.prefix=basename(out.dir))
}

WFUN.quantile <- function(out.dir, WFUN="mean", quant=0.5, deltaG=1, mc.cores=1) {
	data <- WFUN.dir(out.dir, WFUN=WFUN, deltaG=deltaG, mc.cores=mc.cores)
	apply(do.call(rbind, data), 2, quantile, prob=quant)
}

WFUN.quantile.cache <- function(out.dir, WFUN="mean", quant=0.5, deltaG=1, mc.cores=1, cache.subdir=paste0("Rcache-W", WFUN)) {
	cache.fun(WFUN.quantile, out.dir=out.dir, WFUN=WFUN, quant=quant, deltaG=deltaG, mc.cores=mc.cores, cache.subdir=cache.subdir, file.prefix=basename(out.dir))
}
