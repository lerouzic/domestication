#!/usr/bin/env Rscript

# Figure D: evolution of the reaction norm

source("../src/analysis_tools.R")

library(parallel)
mc.cores <- min(12, detectCores()-1)

my.mean.sim <- function(x) mean.sim(x, max.reps=if (detectCores() > 23) Inf else 5, mc.cores=mc.cores) # For tests 
my.mean.norm <- function(x) mean.norm(x, FUN=abs, max.reps=if (detectCores() > 23) Inf else 5, mc.cores=mc.cores) # the FUN=abs is important here! 


# The norm is calculated based on the simulation results using a window:
window.size <- 50 # Number of time points (!= number of generations!)
sliding <- FALSE   # Beware of the computational cost


onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

out.dir.default  <- "../cache/simDefault"
out.dir.nobottle <- "../cache/simNobot"

mean.sim.default  <- my.mean.sim(out.dir.default)
mean.sim.nobottle <- my.mean.sim(out.dir.nobottle)

mean.norm.default  <- my.mean.norm(out.dir.default)
mean.norm.nobottle <- my.mean.norm(out.dir.nobottle)

selpattern.default  <- selectionregime.detect(mean.sim.default)[-1] # The first gene is the environmental signal
selpattern.nobottle <- selectionregime.detect(mean.sim.nobottle)[-1]

Ndyn.default <- get.Ndyn(onerep(out.dir.default))

col <- c(p="red", n="black", s="blue")
lty <- c(default=1, nobottle=2)

pdf("figD.pdf", width=5, height=5)

plot(NULL, xlim=range(mean.sim.default[,"Gen"]), ylim=c(0, 1.2), xlab="Generation", ylab="|reaction norm|")

yy.default.pp <- rowMeans(mean.norm.default[,selpattern.default=="pp"])
yy.default.ps <- rowMeans(mean.norm.default[,selpattern.default=="ps"])
yy.default.pn <- rowMeans(mean.norm.default[,selpattern.default=="pn"])

lines(as.numeric(rownames(mean.norm.default)), yy.default.pp, col=col["p"], lty=lty["default"])
lines(as.numeric(rownames(mean.norm.default)), yy.default.ps, col=col["s"], lty=lty["default"])
lines(as.numeric(rownames(mean.norm.default)), yy.default.pn, col=col["n"], lty=lty["default"])

yy.nobottle.pp <- rowMeans(mean.norm.nobottle[,selpattern.nobottle=="pp"])
yy.nobottle.ps <- rowMeans(mean.norm.nobottle[,selpattern.nobottle=="ps"])
yy.nobottle.pn <- rowMeans(mean.norm.nobottle[,selpattern.nobottle=="pn"])

lines(as.numeric(rownames(mean.norm.nobottle)), yy.nobottle.pp, col=col["p"], lty=lty["nobottle"])
lines(as.numeric(rownames(mean.norm.nobottle)), yy.nobottle.ps, col=col["s"], lty=lty["nobottle"])
lines(as.numeric(rownames(mean.norm.nobottle)), yy.nobottle.pn, col=col["n"], lty=lty["nobottle"])

bottleneck.plot(Ndyn.default, y=0, lwd=2)
selectionchange.plot(mean.sim.default, y=0, cex=1.5)

legend("topright", pch=22, col=col, pt.bg=col, legend=c("Plastic->Pastic", "Plastic->Neutral", "Pastic->Stable"), bty="n", cex=0.8)
legend("topleft", lty=c(1,2), col="darkgray", legend=c("Bottleneck", "No bottleneck"), bty="n", cex=0.8)

dev.off()
