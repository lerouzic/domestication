# Various algorithms to compute phenotypic and genotypic variation


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

molec.variation.neutral <- function(out.table, expr.thresh) {
	allvar <- colnames(out.table)[grep(colnames(out.table), pattern="VarAll")]
	avg.expr <- colMeans(out.table[,grepl(colnames(out.table), pattern="MPhen")])[-1]
	target.genes <- 1 + which(avg.expr <= expr.thresh)
	if(length(target.genes) == 0) return(NULL)
	ans <- t(apply(out.table[,allvar], 1, function(vv) { mv <- matrix(vv, byrow=TRUE, ncol=sqrt(length(vv))); rowMeans(mv[,target.genes,drop=FALSE]) }))
	rownames(ans) <- as.character(out.table[,"Gen"])
	ans
}

molec.variation.neutral.files <- function(out.dir, expr.thresh, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
			tt <- read.table(ff, header=TRUE)
			mvn <- molec.variation.neutral(tt, expr.thresh)
			rm(tt); gc()
			mvn
		}, mc.cores=mc.cores)
	arr <- do.call(abind, c(ans, list(along=3)))
	rowMeans(arr, dims=2)
}

molec.variation.neutral.files.cache <- function(out.dir, expr.thresh, mc.cores=1) {
	cache.fun(molec.variation.neutral.files, out.dir=out.dir, expr.thresh=expr.thresh, mc.cores=mc.cores, cache.subdir="Rcache-vneutral")
}
