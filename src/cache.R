
suppressMessages(library(digest))

default.cache.dir <- "../cache/"

cache.fun <- function(FUN, ..., cache.dir=default.cache.dir, cache.subdir="misc", cache.file=paste0(digest(as.list(...))), ".rds"), rds.version=2) {
	use.cache <- !is.null(cache.dir)
	if (!use.cache) warning ("No cache dir specified. Recalculating...")
	if (use.cache) {
		full.cache.dir <- file.path(cache.dir, cache.subdir)
		if (!dir.exists(full.cache.dir))
			dir.create(full.cache.dir)
		full.cache.file <- file.path(full.cache.dir, cache.file)
		if (file.exists(cache.file))
			return(readRDS(cache.file))
	}
	ans <- FUN(x, ...)
	if (use.cache)
		saveRDS(ans, cache.file, version=rds.version)
	return(ans)
}
