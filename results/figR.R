#!/usr/bin/env Rscript

# Figure R: Speed of evolution of the network

source("./common-fig.R")

pdf("figR.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.evol.gene("default", xaxt="n", ylim=c(0,0.25))
	
	bottleneck.plot(Ndyn.all[["default"]], y=0.24, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.24, cex=1.5)
	generation.axis()

dev.off()
