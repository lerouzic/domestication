# Common functions and calculations for figs B and C

source("../src/analysis_tools.R")

library(parallel)
mc.cores <- min(12, detectCores()-1)

use.cache <- TRUE
cache.dir <- default.cache.dir # defined in ../src/cache.R

window.avg <- 9 # Size of the moving average window
first.gen  <- 0  # First generation for the time series

onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

out.dir.default  <- file.path(cache.dir, "simDefault")
out.dir.nobottle <- file.path(cache.dir, "simNobot")
out.dir.noselc   <- file.path(cache.dir, "simNoselc")
out.dir.nosel    <- file.path(cache.dir, "simNosel")

mean.sim.default  <- mean.sim.cache(out.dir.default)
mean.sim.nobottle <- mean.sim.cache(out.dir.nobottle)
mean.sim.noselc   <- mean.sim.cache(out.dir.noselc)
mean.sim.nosel    <- mean.sim.cache(out.dir.nosel)

selpattern.default  <- selectionregime.detect(mean.sim.default)[-1] # The first gene is the environmental signal
selpattern.nobottle <- selectionregime.detect(mean.sim.nobottle)[-1]
selpattern.noselc   <- selectionregime.detect(mean.sim.noselc)[-1]
selpattern.nosel    <- selectionregime.detect(mean.sim.nosel)[-1]

Ndyn.default    <- get.Ndyn.cache(onerep(out.dir.default))
Ndyn.noselc     <- get.Ndyn.cache(onerep(out.dir.noselc))
Ndyn.nosel      <- get.Ndyn.cache(onerep(out.dir.nosel))

col <- c(c="blue", s="blue", p="red", n="black")
lty <- c(c=1, s=1, p=2, n=3)

