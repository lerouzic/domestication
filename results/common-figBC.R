# Common functions and calculations for figs B and C

source("../src/analysis_tools.R")

library(parallel)
mc.cores <- min(12, detectCores()-1)

mean.sim <- function(out.dir, max.reps=5) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)

	tt <- mclapply(out.files[1:(if (detectCores() > 24 || length(out.files) < max.reps) length(out.files) else max.reps)], read.table, header=TRUE, mc.cores=4) # read.tables on many cores is useless, probably limited by disk speed
	ans <- replicate.mean(tt)
	rm(tt)
	gc()
	return(ans)
}

onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

out.dir.default  <- "../cache/simDefault"
out.dir.nobottle <- "../cache/simNobot"
out.dir.noselc   <- "../cache/simNoselc"
out.dir.nosel    <- "../cache/simNosel"

mean.sim.default  <- mean.sim(out.dir.default)
mean.sim.nobottle <- mean.sim(out.dir.nobottle)
mean.sim.noselc   <- mean.sim(out.dir.noselc)
mean.sim.nosel    <- mean.sim(out.dir.nosel)

selpattern.default  <- selectionregime.detect(mean.sim.default)[-1] # The first gene is the environmental signal
selpattern.nobottle <- selectionregime.detect(mean.sim.nobottle)[-1]
selpattern.noselc   <- selectionregime.detect(mean.sim.noselc)[-1]
selpattern.nosel    <- selectionregime.detect(mean.sim.nosel)[-1]

Ndyn.default    <- get.Ndyn(onerep(out.dir.default))
Ndyn.nobottle   <- get.Ndyn(onerep(out.dir.nobottle))
Ndyn.noselc     <- get.Ndyn(onerep(out.dir.noselc))
Ndyn.nosel      <- get.Ndyn(onerep(out.dir.nosel))

col <- c(c="blue", s="blue", p="red", n="black")
lty <- c(c=1, s=1, p=2, n=3)
