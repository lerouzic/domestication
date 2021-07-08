# Various algorithms to compute phenotypic and genotypic variation

source("../src/analysis_networks.R") # for cleanW

# Phenotypic (expression) variation
pheno.variation <- function(out.table) {
	genes <- colnames(out.table)[grep(colnames(out.table), pattern="VPhen")]
	out.table[,genes]
}

# Molecular (gene network) variation
molec.variation <- function(out.table) {
	allvar <- colnames(out.table)[grep(colnames(out.table), pattern="VarAll")]
	n.genes <- sqrt(length(allvar))
	t(apply(out.table[,allvar], 1, function(x) { rowMeans(matrix(x, ncol=n.genes, byrow=TRUE)) }))
}

molec.variation.dir <- function(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern="VarAll") {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	tt <- results.table(out.files, mc.cores, max.reps, colnames.pattern)
	ans <- lapply(tt, molec.variation)
	rm(tt)
	gc()
	return(ans)
}

molec.variation.mean <- function(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern="VarAll") {
	data <- molec.variation.dir(out.dir, max.reps, mc.cores, colnames.pattern)
	colMeans(do.call(rbind, data))
}

molec.variation.mean.cache <- function(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern="VarAll")
	cache.fun(molec.variation.mean, out.dir=out.dir, max.reps=max.reps, mc.cores=mc.cores, colnames.pattern=colnames.pattern, cache.subdir="Rcache-molvar", file.prefix=basename(out.dir))

molec.variation.quantile <- function(quant=0.5, max.reps=Inf, mc.cores=1, colnames.pattern="VarAll") {
	data <- molec.variation.dir(out.dir, max.reps, mc.cores, colnames.pattern)
	apply(do.call(rbind, data), 2, quantile, prob=quant)
}

molec.variation.quantile.cache <- function(out.dir, quant=0.5, max.reps=Inf, mc.cores=1, colnames.pattern="VarAll")
	cache.fun(molec.variation.quantile, out.dir=out.dir, quand=quant, max.reps=max.reps, mc.cores=mc.cores, colnames.pattern=colnames.pattern, cache.subdir="Rcache-molvar", file.prefix=basename(out.dir))

# Molecular (approximately) neutral variation

molec.variation.neutral.lowexpr <- function(out.table, expr.thresh) {
	allvar <- colnames(out.table)[grep(colnames(out.table), pattern="VarAll")]
	avg.expr <- colMeans(out.table[,grepl(colnames(out.table), pattern="MPhen")])[-1]
	target.genes <- 1 + which(avg.expr <= expr.thresh)
	if(length(target.genes) == 0) return(NULL)
	ans <- t(apply(out.table[,allvar], 1, function(vv) { mv <- matrix(vv, byrow=TRUE, ncol=sqrt(length(vv))); rowMeans(mv[,target.genes,drop=FALSE]) }))
	rownames(ans) <- as.character(out.table[,"Gen"])
	ans
}

model.variation.neutral.cleanW <- function(out.table, expr.thresh, gen=out.table$Gen[nrow(out.table)], env=0.5) {
	print("cleanW")
	allvar <- colnames(out.table)[grep(colnames(out.table), pattern="VarAll")]
	W.table <- out.table[out.table$Gen == gen, grepl(colnames(out.table), pattern="MeanAll")]
	W <- matrix(unlist(W.table, ncol= sqrt(ncol(W.table))), byrow=TRUE)
	cW <- cleanW(W, epsilon=expr.thresh, env=env)
	which.neutral <- which(t(W != cW))            # W is given by row in out.table
	if (length(which.neutral) == 0) return (NULL)
	ans <- out.table[,all.var[which.neutral]]
	rownames(ans) <- as.character(out.table[,"Gen"])
	ans
}

molec.variation.neutral.dir <- function(out.dir, expr.thresh, algorithm=c("lowexpr", "cleanW")[1],  mc.cores=1) {
	stopifnot(algorithm %in% c("lowexpr", "cleanW"))
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
			tt <- read.table(ff, header=TRUE)
			if (algorithm == "lowexpr") {
				mvn <- molec.variation.neutral.lowexpr(tt, expr.thresh)
			} else if (algorithm == "cleanW") {
				mvn <- molec.variation.neutral.cleanW(tt, expr.thresh)
			}
			rm(tt); gc()
			mvn
		}, mc.cores=mc.cores)
	ans
}

# Mean for each gene

molec.variation.neutral.gene.mean <- function(out.dir, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	data <- molec.variation.neutral.dir(out.dir=out.dir, expr.thresh=expr.thresh, mc.cores=mc.cores)
	arr <- do.call(abind, c(data, list(along=3)))
	rowMeans(arr, dims=2)
}

molec.variation.neutral.gene.mean.cache <- function(out.dir, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	cache.fun(molec.variation.neutral.gene.mean, out.dir=out.dir, expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores, cache.subdir="Rcache-vmgneutral", file.prefix=basename(out.dir))
}


# Mean of all genes

molec.variation.neutral.all.mean <- function(out.dir, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	data <- molec.variation.neutral.dir(out.dir=out.dir, expr.thresh=expr.thresh, mc.cores=mc.cores)
	browser()
	arr <- do.call(cbind, lapply(data, function(dd) if (is.null(dd)) NULL else rowMeans(dd)))
	rowMeans(arr)
}

molec.variation.neutral.all.mean.cache <- function(out.dir, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	cache.fun(molec.variation.neutral.all.mean, out.dir=out.dir, expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores, cache.subdir="Rcache-vmaneutral", file.prefix=basename(out.dir))
}

# Quantiles of all genes

molec.variation.neutral.all.quantile <- function(out.dir, quant=0.5, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	data <- molec.variation.neutral.dir(out.dir=out.dir, expr.thresh=expr.thresh, mc.cores=mc.cores)
	arr <- do.call(cbind, lapply(data, function(dd) if (is.null(dd)) NULL else rowMeans(dd)))
	apply(arr, 1, quantile, prob=quant)
}

molec.variation.neutral.all.quantile.cache <- function(out.dir, quant=0.5, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	cache.fun(molec.variation.neutral.all.quantile, out.dir=out.dir, quant=quant, expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores, cache.subdir="Rcache-vqaneutral", file.prefix=basename(out.dir))
}
