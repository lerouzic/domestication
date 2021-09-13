#!/usr/bin/env Rscript

# Figure: evolution of the number of gained/lost connections

source("./common-fig.R")

source("../src/analysis_networks.R")

scenarios <- c("default", "nobot", "noselc")

pdf("fig6A.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.inout.gainloss(scenarios, show.quantiles=TRUE, deltaG=deltaG, ylim=c(-25,25), xaxt="n")

	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=-25, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=-25, cex=1.5)

	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")

	subpanel("A")
dev.off()
