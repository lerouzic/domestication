#!/usr/bin/env Rscript

# Dynamics of the effective pop size

source("./common-fig.R")

mysim <- "default"

pdf("figS1A.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.N(mysim, xaxt="n", xlab="", ylim=c(0,22000))

	generation.axis()
	bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)
	legend("bottomright", lty=1, col=c("black","blue"), legend=c("N", expression(N[e])), cex=cex.legend, bty="n")
	
	subpanel("A")

dev.off()
