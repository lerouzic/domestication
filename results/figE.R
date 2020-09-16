#!/usr/bin/env Rscript

# Figure E: evolution of the network complexity (number of connections)

source("../src/analysis_tools.R")
source("../src/analysis_networks.R")

my.mean.sim <- function(x) mean.sim(x, max.reps=if (detectCores() > 23) Inf else 5, mc.cores=mc.cores) # For tests 
my.mean.connect <- function(x) mean.connect(x, env=0.5, epsilon=NULL, max.reps=if (detectCores() > 23) Inf else 5, mc.cores=mc.cores)


onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

out.dir.default  <- "../cache/simDefault"
out.dir.nobottle <- "../cache/simNobot"

mean.sim.default  <- my.mean.sim(out.dir.default)
mean.sim.nobottle <- my.mean.sim(out.dir.nobottle)

mean.connect.default  <- my.mean.connect(out.dir.default)
mean.connect.nobottle <- my.mean.connect(out.dir.nobottle)

selpattern.default  <- selectionregime.detect(mean.sim.default)[-1] # The first gene is the environmental signal
selpattern.nobottle <- selectionregime.detect(mean.sim.nobottle)[-1]

Ndyn.default <- get.Ndyn(onerep(out.dir.default))
