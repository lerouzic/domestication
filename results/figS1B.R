#!/usr/bin/env Rscript

# Average fitness for default and control simulations

source("./common-fig.R")

scenarios <- c("default", "nobot", "noselc")

pdf("figS1B.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.fitness(scenarios, xaxt="n", xlab="")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=1, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=1, cex=1.5)
	legend("bottomright", lty=lty.sce[scenarios], col="black", legend=legname(scenarios), cex=cex.legend, bty="n")
	subpanel("B")
dev.off()
