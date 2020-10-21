# Common functions and calculations for figs B and C

source("../src/analysis_tools.R")

library(parallel)
mc.cores <- min(12, detectCores()-1)

use.cache <- TRUE
cache.dir <- "../cache"

window.avg <- 9 # Size of the moving average window
first.gen  <- 0  # First generation for the time series

my.mean.sim <- function(x) {
	cache.file <- paste0(file.path(cache.dir, basename(x)), "-mean.rds")
	if (use.cache && file.exists(cache.file)) {
		return(readRDS(cache.file))
	} else {
		ans <- mean.sim(x, max.reps=Inf, mc.cores=mc.cores)
		saveRDS(ans, cache.file, version=2)
		return(ans)
	}
}

onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

out.dir.default  <- "../cache/simDefault"
out.dir.nobottle <- "../cache/simNobot"
out.dir.noselc   <- "../cache/simNoselc"
out.dir.nosel    <- "../cache/simNosel"

mean.sim.default  <- my.mean.sim(out.dir.default)
mean.sim.nobottle <- my.mean.sim(out.dir.nobottle)
mean.sim.noselc   <- my.mean.sim(out.dir.noselc)
mean.sim.nosel    <- my.mean.sim(out.dir.nosel)

selpattern.default  <- selectionregime.detect(mean.sim.default)[-1] # The first gene is the environmental signal
selpattern.nobottle <- selectionregime.detect(mean.sim.nobottle)[-1]
selpattern.noselc   <- selectionregime.detect(mean.sim.noselc)[-1]
selpattern.nosel    <- selectionregime.detect(mean.sim.nosel)[-1]

Ndyn.default    <- get.Ndyn(onerep(out.dir.default))
#~ Ndyn.nobottle   <- get.Ndyn(onerep(out.dir.nobottle))
Ndyn.noselc     <- get.Ndyn(onerep(out.dir.noselc))
Ndyn.nosel      <- get.Ndyn(onerep(out.dir.nosel))

col <- c(c="blue", s="blue", p="red", n="black")
lty <- c(c=1, s=1, p=2, n=3)

