
# Reaction norm (over a given time window)
reaction.norm <- function(partial.out.table, signal.gene="MPhen1") {
	.reac.norm <- function(x, signal) coef(lm(x ~ signal))[2]
	expr <- partial.out.table[,grepl(pattern="MPhen", colnames(partial.out.table))]
	sig <- expr[,signal.gene]
	apply(expr[,colnames(expr) != signal.gene], 2, .reac.norm, signal=sig)
}

# Reaction norm dynamics (for a given window size)
reaction.norm.dyn <- function(out.table, window.size=floor(nrow(out.table)/100), sliding=TRUE, signal.gene="MPhen1") {
	stopifnot(window.size > 1, window.size <= nrow(out.table))
	list.gens <- if (sliding) {
					lapply(1:(nrow(out.table)-window.size+1), function(i.start) i.start:(i.start+window.size-1))
				} else {
					lapply(1:floor(nrow(out.table)/window.size), function(i.start) (1+(i.start-1)*window.size):(i.start*window.size))
				}
	ans <- do.call(rbind, lapply(list.gens, function(x) reaction.norm(out.table[x,], signal.gene=signal.gene)))
	rownames(ans) <- as.character(sapply(list.gens, function(x) round(mean(out.table[x,"Gen"]))))
	ans
}

# Average out all reaction norms from a directory (use FUN=abs to get the absolute value of the norm)
reaction.norm.dir <- function(out.dir, max.reps=Inf, FUN.to.apply=identity, mc.cores=1, sliding=TRUE, window.size=10) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	tt <- results.table(out.files, mc.cores, max.reps)
	nn <- mclapply(tt, reaction.norm.dyn, window.size=window.size, sliding=sliding, mc.cores=mc.cores)
	ans <- lapply(nn, FUN.to.apply)
	rm(tt)
	gc()
	return(ans)
}

reaction.norm.mean <- function(out.dir, max.reps=Inf, FUN.to.apply=identity, mc.cores=1, sliding=TRUE, window.size=10) {
	data <- reaction.norm.dir(out.dir, max.reps=max.reps, FUN.to.apply=FUN.to.apply, mc.cores=mc.cores, sliding=sliding, window.size=window.size)
	replicate.mean(data)
}

reaction.norm.mean.cache <- function(out.dir, max.reps=Inf, FUN.to.apply=identity, mc.cores=1, sliding=TRUE, window.size=10) {
	cache.fun(mean.norm, out.dir=out.dir, max.reps=max.reps, FUN.to.apply=FUN.to.apply, mc.cores=mc.cores, sliding=sliding, window.size=window.size, cache.subdir="Rcache-norm")
}

reaction.norm.quantile <- function(out.dir, quant=0.5, max.reps=Inf, FUN.to.apply=identity, mc.cores=1, sliding=TRUE, window.size=10) {
	data <- reaction.norm.dir(out.dir, max.reps=max.reps, FUN.to.apply=FUN.to.apply, mc.cores=mc.cores, sliding=sliding, window.size=window.size)
	replicate.quantile(data, quant=quant)
}

reaction.norm.quantile.cache <- function(out.dir, quant=0.5, max.reps=Inf, FUN.to.apply=identity, mc.cores=1, sliding=TRUE, window.size=10) {
	cache.fun(reaction.norm.quantile, out.dir=out.dir, quant=quant, max.reps=max.reps, FUN.to.apply=FUN.to.apply, mc.cores=mc.cores, sliding=sliding, window.size=window.size, cache.subdir="Rcache-norm")
}
