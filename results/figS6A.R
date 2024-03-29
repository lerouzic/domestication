#!/usr/bin/env Rscript

# Figure: Evolution of genetic variance-covariance

source("./common-fig.R")
source("../src/analysis_networks.R")

scenarios <- c("default","nobot","noselc")

pdf("figS6A.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)

	plot.Gdiff(scenarios, deltaG=deltaG, xaxt="n", ylim=c(0,0.6))
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
	
	legend(x="topright",legend = legname(scenarios), lty=lty.sce[scenarios], col=col.sce[scenarios], bty="n", cex=cex.legend)
	subpanel("A")
dev.off()
