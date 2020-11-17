#!/usr/bin/env Rscript

# Figure L: Evolution of genetic variance-covariance

library(ade4) # for mantel.rtest

source("./common-fig.R")

source("../src/analysis_networks.R")


scenarios <- c("default","nobot","noselc")

pdf("figL.pdf", width=8, height=4)
	layout(rbind(1:2))
	
	plot.Gdiff(scenarios, deltaG=deltaG, xaxt="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
	
	legend(x="topright",legend = legname(scenarios), lty=lty.sce[scenarios], col=col.sce[scenarios], bty="n")
	subpanel("A")
	
	
	plot.GPC(scenarios, PC=1, xaxt="n")	
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)

	subpanel("B")

dev.off()



