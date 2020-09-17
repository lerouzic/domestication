#!/usr/bin/env Rscript

# Figure E: evolution of the network complexity (number of connections)

# Computing the number of connections for each time point is time consuming. 
# However, a cache system (in ../cache/distW) makes further calls faster.

library(parallel)
mc.cores <- min(12, detectCores()-1)

source("../src/analysis_tools.R")
source("../src/analysis_networks.R")

connect.threshold <- 0.1

my.mean.sim <- function(x) mean.sim(x, max.reps=if (detectCores() > 23) Inf else 5, mc.cores=mc.cores) # For tests 
my.mean.connect <- function(x) mean.connect(x, env=0.5, epsilon=connect.threshold, max.reps=if (detectCores() > 23) Inf else 5, mc.cores=mc.cores)


onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

out.dir.default  <- "../cache/simDefault"
out.dir.nobottle <- "../cache/simNobot"
out.dir.noselc   <- "../cache/simNoselc"

mean.sim.default  <- my.mean.sim(out.dir.default)
#~ mean.sim.nobottle <- my.mean.sim(out.dir.nobottle)
#~ mean.sim.noselc   <- my.mean.sim(out.dir.noselc)

mean.connect.default  <- my.mean.connect(out.dir.default)
mean.connect.nobottle <- my.mean.connect(out.dir.nobottle)
mean.connect.noselc   <- my.mean.connect(out.dir.noselc)

selpattern.default  <- selectionregime.detect(mean.sim.default)[-1] # The first gene is the environmental signal
#~ selpattern.nobottle <- selectionregime.detect(mean.sim.nobottle)[-1]
#~ selpattern.noselc   <- selectionregime.detect(mean.sim.noselc)[-1]

Ndyn.default <- get.Ndyn(onerep(out.dir.default))

col <- c(default="black", nobottle="darkgreen", noselc="orange")
leg <- c(default="Bottleneck + sel change", nobottle="Sel change (no bottleneck)", noselc="Bottleneck (no sel change)")

pdf("figE.pdf", width=5, height=5)

plot(NULL, xlim=range(mean.sim.default[,"Gen"]), ylim=c(0, max(c(mean.connect.default, mean.connect.nobottle, mean.connect.noselc))), xlab="Generation", ylab="Nb connections")

lines(as.numeric(names(mean.connect.default)),  mean.connect.default,  col=col["default"])
lines(as.numeric(names(mean.connect.nobottle)), mean.connect.nobottle, col=col["nobottle"])
lines(as.numeric(names(mean.connect.noselc)),   mean.connect.noselc,   col=col["noselc"])

bottleneck.plot(Ndyn.default, y=0, lwd=2)
selectionchange.plot(mean.sim.default, y=0, cex=1.5)

legend("topleft", pch=22, col=col, pt.bg=col, legend=leg, bty="n", cex=0.8)

dev.off()
