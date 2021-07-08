#!/usr/bin/env Rscript

# Figure: evolution of the number of gained/lost connections

source("./common-fig.R")

source("../src/analysis_networks.R")

scenarios <- c("default", "nobot", "noselc")

pdf("fig5A.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.inout.gainloss(scenarios, deltaG=deltaG, ylim=c(-20, 20), xaxt="n")
	
	legend("topright", lty=c(1,1,lty.sce[scenarios]), col=c(col.gl, rep("black", length(scenarios))), legend=legname(c(names(col.gl), scenarios)), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=-19.5, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=-19.5, cex=1.5)

	subpanel("A")
dev.off()
