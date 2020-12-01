#!/usr/bin/env Rscript

# Figure S: Molecular variance of neutral sites

source("./common-fig.R")

pdf("figS.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	ylab <- "Molecular variance"
	y.factor <- c('(""%*% 10^{-4})' = 10000)
	if (y.factor != 1) ylab <- parse(text=paste0('"', ylab, ' "*', names(y.factor)))
	ylim <- c(0, 2e-4)*y.factor
	
	plot.var.neutral.gene("default", xaxt="n", ylab=ylab, ylim=ylim, y.factor=y.factor)
	
	bottleneck.plot(Ndyn.all[["default"]], y=0.24, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.24, cex=1.5)
	generation.axis()

dev.off()
