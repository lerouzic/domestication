# Common functions and calculations for figs B and C

source("../src/analysis_tools.R")

library(parallel)
mc.cores <- min(12, detectCores()-1)

cache.dir <- default.cache.dir # defined in ../src/cache.R

window.avg <- 9 # Size of the moving average window
first.gen  <- 0  # First generation for the time series

onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

out.dir.default  <- file.path(cache.dir, "simDefault")
out.dir.nobottle <- file.path(cache.dir, "simNobot")
out.dir.noselc   <- file.path(cache.dir, "simNoselc")
out.dir.nosel    <- file.path(cache.dir, "simNosel")

tokeep <- "Mphen|FitOpt|MFit|VFit"

mean.sim.default  <- mean.sim.cache(out.dir.default, colnames.pattern=tokeep)
mean.sim.nobottle <- mean.sim.cache(out.dir.nobottle, colnames.pattern=tokeep)
mean.sim.noselc   <- mean.sim.cache(out.dir.noselc, colnames.pattern=tokeep)
mean.sim.nosel    <- mean.sim.cache(out.dir.nosel, colnames.pattern=tokeep)

selpattern.default  <- selectionregime.detect(mean.sim.default)[-1] # The first gene is the environmental signal
selpattern.nobottle <- selectionregime.detect(mean.sim.nobottle)[-1]
selpattern.noselc   <- selectionregime.detect(mean.sim.noselc)[-1]
selpattern.nosel    <- selectionregime.detect(mean.sim.nosel)[-1]

Ndyn.default    <- get.Ndyn.cache(onerep(out.dir.default))
Ndyn.noselc     <- get.Ndyn.cache(onerep(out.dir.noselc))
Ndyn.nosel      <- get.Ndyn.cache(onerep(out.dir.nosel))

col <- c(c="blue", s="blue", p="red", n="black")
lty <- c(c=1, s=1, p=2, n=3)

