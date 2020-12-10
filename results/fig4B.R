#!/usr/bin/env Rscript

# Figure: Evolution of average genetic correlations 

source("./common-fig.R")

source("../src/analysis_networks.R")

scenarios <- c("default", "nobot", "noselc")

pdf("fig4B.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.Gcor(scenarios, ylim=c(0,0.2), xaxt="n")
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0.05, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.05, cex=1.5)

dev.off()
