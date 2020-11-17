#!/usr/bin/env Rscript

# Figure G: number of modules

source("./common-fig.R")

connect.threshold <- 0.1
env <- 0.5
directed <- FALSE
use.cache <- TRUE

source("../src/analysis_networks.R")

scenarios <- c("default","nobot","noselc")

pdf("figG.pdf", width=10, height=5)
	
	layout(t(1:2))
	
	plot.network.feature(scenarios, what="nbconn", xaxt="n", xlab="")

	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
	legend("topleft", lty=lty.sce[scenarios], col="black", legend=legname(scenarios), bty="n", cex=0.8)
	legend("topright", pch=17, col=col.algo, legend=names(col.algo), bty="n", cex=0.8)
	
	subpanel("A")


	plot.network.feature(scenarios, what="modularity", xaxt="n", xlab="", ylim=c(0, 0.5))

	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
	legend("topleft", lty=lty.sce[scenarios], col="black", legend=legname(scenarios), bty="n", cex=0.8)
	legend("topright", pch=17, col=col.algo, legend=names(col.algo), bty="n", cex=0.8)
	
	subpanel("B")

dev.off()
