#!/usr/bin/env Rscript

# Figure E: evolution of the network complexity (number of connections)

# Computing the number of connections for each time point is time consuming. 
# However, a cache system (in ../cache/distW) makes further calls faster.

source("./commonfig.R")

source("../src/analysis_networks.R")

connect.threshold <- 0.1

mean.connect.default  <- mean.connect.cache(out.dir.default)
mean.connect.nobottle <- mean.connect.cache(out.dir.nobottle)
mean.connect.noselc   <- mean.connect.cache(out.dir.noselc)

my.mean.connect.default <-  mov.avg(mean.connect.default,  as.numeric(names(mean.connect.default)),  size=window.avg, min.gen=first.gen)
my.mean.connect.nobottle <- mov.avg(mean.connect.nobottle, as.numeric(names(mean.connect.nobottle)), size=window.avg, min.gen=first.gen)
my.mean.connect.noselc <-   mov.avg(mean.connect.noselc,   as.numeric(names(mean.connect.noselc)),   size=window.avg, min.gen=first.gen)

scenarios <- c("default", "nobottle", "noselc")

pdf("figE.pdf", width=5, height=5)

plot(NULL, 
	xlim=c(first.gen, max(mean.sim.default[,"Gen"])), 
	ylim=c(0, max(c(mean.connect.default, mean.connect.nobottle, mean.connect.noselc))), 
	xlab="Generation", ylab="Nb connections", xaxt="n")

lines(as.numeric(names(my.mean.connect.default)),  my.mean.connect.default,  col=col.sce["default"],  lty=lty.sce["default"])
lines(as.numeric(names(my.mean.connect.nobottle)), my.mean.connect.nobottle, col=col.sce["nobottle"], lty=lty.sce["nobottle"])
lines(as.numeric(names(my.mean.connect.noselc)),   my.mean.connect.noselc,   col=col.sce["noselc"],   lty=lty.sce["noselc"])

bottleneck.plot(Ndyn.default, y=0, lwd=2)
selectionchange.plot(mean.sim.default, y=0, cex=1.5)

legend("topleft", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), bty="n", cex=0.8)
generation.axis()

dev.off()
