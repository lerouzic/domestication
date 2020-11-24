#!/usr/bin/env Rscript

# Figure L: Evolution of genetic variance-covariance

source("./common-fig.R")

source("../src/analysis_networks.R")

scenarios <- c("default","nobot","noselc")

pdf("figL.pdf", width=2*panel.width, height=panel.height)
	layout(rbind(1:2))

	plot.GPC(scenarios, PC=1, xaxt="n")	
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)

	subpanel("A")
	
	
	plot.Grank(scenarios, xaxt="n")	
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)

	subpanel("B")
	
dev.off()



