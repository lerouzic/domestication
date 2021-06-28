# Number of network connections

number.connections <- function(W, ...) {
	cW <- cleanW.cache(W=W, ...)
	sum(cW != 0)
}

number.connections.dyn <- function(out.table) { 
	W.table <- out.table[,grepl(colnames(out.table), pattern="MeanAll")]
	
	net.size <- sqrt(ncol(W.table))
	nb.conn <- sapply(1:nrow(W.table), function(i) { 
			W <- matrix(unlist(W.table[i,]), ncol=net.size, byrow=TRUE)
			if(nrow(W) != ncol(W)) return(NA)
			ans <- try(number.connections(W))
			if (class(ans) == "try-error") NA else ans
		})
	names(nb.conn) <- as.character(out.table[,"Gen"])
	nb.conn
}

number.connections.dir <- function(out.dir, max.reps=Inf, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	tt <- results.table(out.files, mc.cores, max.reps)
	ans <- mclapply(tt, number.connections.dyn, mc.cores=mc.cores)
	rm(tt)
	gc()
	return(ans)
}

# Average out all network connections from a directory 
number.connections.mean <- function(out.dir, max.reps=Inf, mc.cores=1) {
	data <- number.connections.dir(out.dir, max.reps, mc.cores=mc.cores)
	colMeans(do.call(rbind, nn), na.rm=TRUE)
}

number.connections.mean.cache <- function(out.dir, max.reps=Inf, mc.cores=1) {
	cache.fun(number.connections.mean, out.dir=out.dir, max.reps=max.reps, mc.cores=mc.cores, cache.subdir="Rcache-connect")
}

# Quantiles
number.connections.quantile <- function(out.dir, quant=0.5, max.reps=Inf, mc.cores=1) {
	data <- number.connections.dir(out.dir, max.reps, mc.cores=mc.cores)
	apply(do.call(rbind, nn), 2, quantile, prob=quant)
}

number.connections.quantile.cache <- function(out.dir, quant=0.5, max.reps=Inf, mc.cores=1) {
	cache.fun(number.connections.quantile, quant=quant, out.dir=out.dir, max.reps=max.reps, mc.cores=mc.cores, cache.subdir="Rcache-connect")
}


################################################# In/Out connections

inout.connections <- function(W, ...) {
	cW <- cleanW.cache(W=W, ...)
	list(connect.in=rowSums(cW != 0), connect.out=colSums(cW != 0))
}

# in/out connections at a specific generation
inout.gen <- function(out.dir, gen, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		if (!gen %in% tt$Gen) return(NA)
		W <- tt[tt[,"Gen"] == gen, grepl(colnames(tt), pattern="MeanAll")]
		rm(tt); gc()
		W <- matrix(unlist(W), ncol=sqrt(length(W)), byrow=TRUE)
		inout.connections(W)
	}, mc.cores=mc.cores)
	ans[!sapply(ans, function(x) length(x) < 2 || is.na(x))]
}

inout.gen.cache <- function(out.dir, gen, mc.cores=1) {
	cache.fun(inout.gen, out.dir=out.dir, gen=gen, mc.cores=mc.cores, cache.subdir="Rcache-inout")
}

delta.inout <- function(W, Wref) {
	cW <- cleanW.cache(W)
	cWref <- cleanW.cache(Wref)
	
	c(gain=sum(cWref == 0 & cW != 0), loss=sum(cWref != 0 & cW == 0))
}

# Returns two columns: one with the gains, another with the losses
delta.inout.dyn <- function(out.table, deltaG=1, mc.cores=1) {
	gen <- out.table[,"Gen"]
	seqgen <- seq(1, length(gen), length.out=min(length(gen), 1+gen[length(gen)] %/% deltaG))
	listW <- Wlist.table(out.table[seqgen,])
	
	ans <- mclapply(1:(length(listW)-1), function(i) delta.inout(listW[[i+1]], listW[[i]]), mc.cores=mc.cores)
	ans <- do.call(rbind, ans)
	rownames(ans) <- as.character(gen[seqgen])[-1]
	ans
}

delta.inout.dir <- function(out.dir, deltaG=NA, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		delta.inout.dyn(tt, deltaG, mc.cores=1)
	}, mc.cores=mc.cores)
	l.ans <- sapply(ans, nrow)
	ans <- ans[l.ans == max(l.ans)] # removing incomplete data sets (ongoing simulations)
	ans
}

delta.inout.mean <- function(out.dir, deltaG=NA, mc.cores=1) {
	data <- delta.inout.dir(out.dir, deltaG=deltaG, mc.cores=mc.cores)
	replicate.mean(data)
}

delta.inout.mean.cache <- function(out.dir, deltaG=NA, mc.cores=1) {
	cache.fun(delta.inout.mean, out.dir=out.dir, deltaG=deltaG, mc.cores=mc.cores, cache.subdir="Rcache-inout")
}

delta.inout.quantile <- function(out.dir, quant=0.5, deltaG=NA, mc.cores=1) {
	data <- delta.inout.dir(out.dir, deltaG=deltaG, mc.cores=mc.cores)
	replicate.quantile(data, quant=quant)
}

delta.inout.quantile.cache <- function(out.dir, quant=0.5, deltaG=NA, mc.cores=1) {
	cache.fun(delta.inout.quantile, out.dir=out.dir, quant=quant, deltaG=deltaG, mc.cores=mc.cores, cache.subdir="Rcache-inout")
}


