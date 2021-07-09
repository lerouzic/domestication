# Various algorithms to compute phenotypic and genotypic variation

source("../src/analysis_tools.R")
source("../src/analysis_networks.R") # for cleanW

# Phenotypic (expression) variation
pheno.variation <- function(out.table) {
	genes <- colnames(out.table)[grep(colnames(out.table), pattern="VPhen")]
	out.table[,genes]
}

pheno.variation.dir <- function(out.dir,  max.reps=Inf, mc.cores=1) {
	out.files <- out.files(out.dir)
	tt <- results.table(out.files, mc.cores, max.reps, colnames.pattern="VPhen")
	ans <- mclapply(tt, pheno.variation, mc.cores=mc.cores)
	rm(tt)
	gc()
	return(ans)
}

# Mean phenotypic (expression) variance over all genes (except the first "env"). 
pheno.variation.mean <- function(out.dir, max.reps=Inf, mc.cores=1, without="VPhen1") {
	data <- pheno.variation.dir(out.dir, max.reps, mc.cores)
	ans <- colMeans(do.call(rbind, lapply(data, function(x) rowMeans(x[,!(colnames(x) %in% without)]))))
	names(ans) <- rownames(data[[1]])
	ans
}

pheno.variation.mean.cache <- function(out.dir, max.reps=Inf, mc.cores=1)
	cache.fun(pheno.variation.mean, out.dir=out.dir, max.reps=max.reps, mc.cores=mc.cores, cache.subdir="Rcache-phenvar", file.prefix=basename(out.dir))
	
pheno.variation.quantile <- function(out.dir, quant=0.5, max.reps=Inf, mc.cores=1) {
	data <- pheno.variation.dir(out.dir, max.reps, mc.cores)
	apply(do.call(rbind, lapply(data, rowMeans)), 2, quantile, prob=quant)
}

pheno.variation.quantile.cache <- function(out.dir, quant=0.5, max.reps=Inf, mc.cores=1, colnames.pattern="VPhen")
	cache.fun(pheno.variation.quantile, out.dir=out.dir, quant=quant, max.reps=max.reps, mc.cores=mc.cores, colnames.pattern=colnames.pattern, cache.subdir="Rcache-phenvar", file.prefix=basename(out.dir))


# Molecular (gene network) variation
molec.variation <- function(out.table) {
	allvar <- colnames(out.table)[grep(colnames(out.table), pattern="VarAll")]
	n.genes <- sqrt(length(allvar))
	t(apply(out.table[,allvar], 1, function(x) { rowMeans(matrix(x, ncol=n.genes, byrow=TRUE)) }))
}

molec.variation.dir <- function(out.dir, max.reps=Inf, mc.cores=1, colnames.pattern="VarAll") {
	out.files <- out.files(out.dir)
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
	which.neutral.genes <- 1 + which(avg.expr <= expr.thresh)
	mm <- matrix(FALSE, ncol=sqrt(length(allvar)), nrow=sqrt(length(allvar)))
	mm[,which.neutral.genes] <- TRUE
	which.neutral.W <- which(c(mm))
	ans <- out.table[,allvar]
	ans[,-which.neutral.W] <- NA
	rownames(ans) <- as.character(out.table[,"Gen"])
	ans
}

molec.variation.neutral.cleanW <- function(out.table, expr.thresh, gen=out.table$Gen[nrow(out.table)], env=0.5) {
	allvar <- colnames(out.table)[grep(colnames(out.table), pattern="VarAll")]
	W.table <- out.table[out.table$Gen == gen, grepl(colnames(out.table), pattern="MeanAll")]
	W <- matrix(unlist(W.table), ncol= sqrt(ncol(W.table)), byrow=TRUE)
	cW <- cleanW(W, epsilon=expr.thresh, env=env)
	which.neutral.W <- which(t(W != cW))            # W is given by row in out.table
	ans <- out.table[,allvar]
	ans[,-which.neutral.W] <- NA
	rownames(ans) <- as.character(out.table[,"Gen"])
	ans
}

molec.variation.neutral.dir <- function(out.dir, expr.thresh, algorithm=c("lowexpr", "cleanW")[1],  mc.cores=1) {
	stopifnot(algorithm %in% c("lowexpr", "cleanW"))
	out.files <- out.files(out.dir)
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
	data.alleles <- molec.variation.neutral.dir(out.dir=out.dir, expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores)
	data.genes <- lapply(data.alleles, function(dd) t(apply(dd, 1, function(aa) {mm <- matrix(unlist(aa), ncol=sqrt(length(aa))); rowMeans(mm, na.rm=TRUE)})))
	arr <- do.call(abind, c(data.genes, list(along=3)))
	rowMeans(arr, dims=2, na.rm=TRUE)
}

molec.variation.neutral.gene.mean.cache <- function(out.dir, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	cache.fun(molec.variation.neutral.gene.mean, out.dir=out.dir, expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores, cache.subdir="Rcache-vmgneutral", file.prefix=basename(out.dir))
}


# Mean of all genes

molec.variation.neutral.all.mean <- function(out.dir, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	data <- molec.variation.neutral.dir(out.dir=out.dir, expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores)
	arr <- do.call(rbind, lapply(data, function(dd) if (is.null(dd)) NULL else rowMeans(dd, na.rm=TRUE)))
	ans <- colMeans(arr, na.rm=TRUE)
	names(ans) <- rownames(data[[1]])
	ans
}

molec.variation.neutral.all.mean.cache <- function(out.dir, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	cache.fun(molec.variation.neutral.all.mean, out.dir=out.dir, expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores, cache.subdir="Rcache-vmaneutral", file.prefix=basename(out.dir))
}

# Quantiles of all genes

molec.variation.neutral.all.quantile <- function(out.dir, quant=0.5, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	data <- molec.variation.neutral.dir(out.dir=out.dir, expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores)
	arr <- do.call(cbind, lapply(data, function(dd) if (is.null(dd)) NULL else rowMeans(dd, na.rm=TRUE)))
	ans <- apply(arr, 1, quantile, prob=quant, na.rm=TRUE)
	names(ans) <- rownames(data[[1]])
	ans
}

molec.variation.neutral.all.quantile.cache <- function(out.dir, quant=0.5, expr.thresh, algorithm=c("lowexpr", "cleanW")[1], mc.cores=1) {
	cache.fun(molec.variation.neutral.all.quantile, out.dir=out.dir, quant=quant, expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores, cache.subdir="Rcache-vqaneutral", file.prefix=basename(out.dir))
}
