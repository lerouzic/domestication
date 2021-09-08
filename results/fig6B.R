#!/usr/bin/env Rscript

# Figure: Evolution of the number of clusters

source("./common-fig.R")

source("../src/analysis_networks.R")

scenarios <- c("default", "nobot", "noselc")
dG <- 1000

pdf("fig6B.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.network.feature(scenarios, what="clusters", xaxt="n", deltaG=dG, ylab="Number of clusters", ylim=c(6, 12))
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=11.5, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=11.5, cex=1.5)
	subpanel("B")
dev.off()
