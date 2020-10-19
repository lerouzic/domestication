#!/usr/bin/env Rscript

# First figure: dynamics of the effective pop size and average fitness +/- sd for a random population

library(parallel)
mc.cores <- min(12, detectCores()-1)

# The analysis is restricted to 5 files when not on the server
my.mean.sim <- function(x) mean.sim(x, max.reps=if (detectCores() > 23) Inf else 5)
my.var.sim  <- function(x) var.sim (x, max.reps=if (detectCores() > 23) Inf else 5)

onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

source("../src/analysis_tools.R")

out.dir <- "../cache/simDefault"

mean.sim  <- my.mean.sim(out.dir)
var.sim   <- my.var.sim(out.dir)

mfit <- mean.sim[,"MFit"]
vfit <- mean.sim[,"VFit"]
gen <-  mean.sim[,"Gen"]

N <- get.Ndyn(onerep(out.dir))[gen+1] 

pdf("figA.pdf", width=8, height=4)

layout(t(1:2))

plot(gen, N, type="l", xlim=range(gen), ylim=c(0, max(N)), ylab="Population size", xlab="Generations")
lines(gen, N/(1+4*vfit/mfit/mfit), col="blue")
bottleneck.plot(N, y=1, lwd=2)
selectionchange.plot(mean.sim, y=1, cex=1.5)
legend("bottomright", lty=1, col=c("black","blue"), legend=c("N", expression(N[e])), bty="n")

plot(NULL, xlab="Generations", ylab="Fitness", xlim=range(gen), ylim=c(0,1))
#~ lines(gen, mfit+sqrt(vfit), lty=1, col="blue")
#~ lines(gen, mfit-sqrt(vfit), lty=1, col="blue")
lines(gen, mfit, lwd=2)
bottleneck.plot(N, y=1, lwd=2)
selectionchange.plot(mean.sim, y=1, cex=1.5)

dev.off()
