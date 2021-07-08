#!/usr/bin/env Rscript

# Figure: Speed of evolution of the network

source("./common-fig.R")

scenarios <- c("default", "nobot", "noselc")

pdf("fig2A.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	sel.pat <- substr(selpattern.all[["default"]], 1, 1)
	sel.pat[sel.pat == "c"] <- "s" # No need to distinguish constant and stable
	
	plot.evol(scenarios, xaxt="n", ylim=c(0,0.15))
	
	bottleneck.plot(Ndyn.all[["default"]], y=0.24, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.24, cex=1.5)
	generation.axis()
	
	legend("topright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), bty="n", cex=cex.legend)

	subpanel("A")
dev.off()
