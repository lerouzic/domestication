#!/usr/bin/env Rscript

# Figure R: Speed of evolution of the network

source("./common-fig.R")

pdf("figR.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.evol.gene("default", xaxt="n")
	
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
	generation.axis()

dev.off()
