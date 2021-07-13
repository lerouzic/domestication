
######## Evolution of the G matrix

Gdiff <- function(G, G.ref) {
	# A bit of cleaning is necessary : first generation and "environmental" gene can mess up the vcov
	if (any(!is.finite(G)) || any(!is.finite(G.ref))) return(NA)
	diag(G)[diag(G) < 1e-8] <- 1e-8
	diag(G.ref)[diag(G.ref) < 1e-8] <- 1e-8
	# covtransf is defined in common-precalc.R, it turns vcov into distance matrices
	mt <- try(mantel.rtest(covtransf(G), covtransf(G.ref), nrepet=1)$obs)
	if (class(mt) == "try-error") mt <- 1
	1 - mt
}

Gdiff.dyn <- function(out.table, deltaG=NA, mc.cores=1) {
	gen <- out.table[,"Gen"]
	seqgen <- seq(1, length(gen), length.out=min(length(gen), 1+gen[length(gen)] %/% deltaG))
	listG <- Glist.table(out.table[seqgen,])
	
	ans <- mclapply(1:(length(listG)-1), function(i) Gdiff(listG[[i+1]], listG[[i]]), mc.cores=mc.cores)
	ans <- do.call(c, ans)
	names(ans) <- as.character(gen[seqgen])[-1]
	ans	
}

Gdiff.dir <- function(out.dir, deltaG=NA, mc.cores=1) {
	out.files <- out.files(out.dir)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		Gdiff.dyn(tt, deltaG, mc.cores=1)
	}, mc.cores=mc.cores)

	ansl <- sapply(ans, length)
	ans <- ans[ansl == max(ansl)]
	ans
}

Gdiff.mean <- function(out.dir, deltaG=NA, mc.cores=1) {
	data <- Gdiff.dir(out.dir, deltaG=deltaG, mc.cores=mc.cores)
	colMeans(do.call(rbind, data))
}

Gdiff.mean.cache <- function(out.dir, deltaG=NA, mc.cores=1) {
	cache.fun(Gdiff.mean, out.dir=out.dir, deltaG=deltaG, mc.cores=mc.cores, cache.subdir="Rcache-Gdiff", file.prefix=basename(out.dir))
}

Gdiff.quantile <- function(out.dir, quant=0.5, deltaG=NA, mc.cores=1) {
	data <- Gdiff.dir(out.dir, deltaG=deltaG, mc.cores=mc.cores)
	apply(do.call(rbind, data), 2, quantile, prob=quant)
}

Gdiff.quantile.cache <- function(out.dir, quant=0.5, deltaG=NA, mc.cores=1) {
	cache.fun(Gdiff.quantile, out.dir=out.dir, quant=quant, deltaG=deltaG, mc.cores=mc.cores, cache.subdir="Rcache-Gdiff", file.prefix=basename(out.dir))
}

# Evolution of the average genetic correlation

Gcor.dyn <- function(out.table) {
	.meancor <- function(G) {
		diag(G)[diag(G) < 1e-12] <- 1e-12 # Prevent numerical errors 
		C <- cov2cor(G)
		diag(C) <- NA # We don't want to count the diagonal when computing the mean
		mean(abs(C), na.rm=TRUE)
	}
	Glist <- Glist.table(out.table)
	sapply(Glist, .meancor)
}


Gcor.dir <- function(out.dir, mc.cores=1) {
	out.files <- out.files(out.dir)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		gc <- Gcor.dyn(tt)
		rm(tt)
		gc
	}, mc.cores=mc.cores)
	
	ansl <- sapply(ans, length)
	ans <- ans[ansl==max(ansl)]
	ans
}

Gcor.mean <- function(out.dir, mc.cores=1) {
	data <- Gcor.dir(out.dir, mc.cores)
	colMeans(do.call(rbind, data), na.rm=TRUE)
}

Gcor.mean.cache <- function(out.dir, mc.cores=1) {
	cache.fun(Gcor.mean, out.dir=out.dir, mc.cores=mc.cores, cache.subdir="Rcache-Gcor", file.prefix=basename(out.dir))
}

Gcor.quantile <- function(out.dir, quant=0.5, mc.cores=1) {
	data <- Gcor.dir(out.dir, mc.cores)
	apply(do.call(rbind, data), 2, quantile, prob=quant, na.rm=TRUE)
}

Gcor.quantile.cache <- function(out.dir, quant=0.5, mc.cores=1) {
	cache.fun(Gcor.quantile, out.dir=out.dir, quant=quant, mc.cores=mc.cores, cache.subdir="Rcache-Gcor", file.prefix=basename(out.dir))
}
