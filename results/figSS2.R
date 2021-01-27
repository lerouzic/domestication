#!/usr/bin/env Rscript

# Figure: Evolution of various graph features

source("./common-fig.R")

source("../src/analysis_networks.R")

scenarios <- c("default", "nobot", "noselc")
dG <- 1

pdf("figSS2a.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.network.feature(scenarios, what="path", xaxt="n", deltaG=dG)
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0.05, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.05, cex=1.5)
dev.off()

pdf("figSS2b.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.network.feature(scenarios, what="diameter", xaxt="n", deltaG=dG)
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0.05, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.05, cex=1.5)
dev.off()

pdf("figSS2c.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.network.feature(scenarios, what="clusters", xaxt="n", deltaG=dG)
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0.05, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.05, cex=1.5)
dev.off()
