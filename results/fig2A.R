#!/usr/bin/env Rscript

# Figure: Speed of evolution of the network

source("./common-fig.R")

pdf("fig2A.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	sel.pat <- substr(selpattern.all[["default"]], 1, 1)
	sel.pat[sel.pat == "c"] <- "s" # No need to distinguish constant and stable
	
	plot.evol.gene("default", xaxt="n", ylim=c(0,0.25))
	
	bottleneck.plot(Ndyn.all[["default"]], y=0.24, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.24, cex=1.5)
	generation.axis()
	
	legend("topright", lty=c(0, 1, 1, 1), col=c(0, col.sel[unique(sel.pat)]), legend=c("Before dom:", "Stable","Plastic", "Non-selected"), cex=cex.legend, bty="n")
	legend("right", lty=c(0, lty.sel[unique(sel.pat)]), col=c(0, 1, 1, 1), legend=c("After dom:", "Stable","Plastic", "Non-selected"), cex=cex.legend, bty="n")


	subpanel("A")
dev.off()
