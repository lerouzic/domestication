#!/usr/bin/env Rscript

# Figure S: Molecular variance of neutral sites

source("./common-fig.R")

pdf("figS.pdf", width=3*panel.width, height=panel.height)
	layout(t(1:3))
	par(mar=mar.notitle)
	
	ylab <- "Molecular variance"
	y.factor <- c('(""%*% 10^{-4})' = 10000)
	if (y.factor != 1) ylab <- parse(text=paste0('"', ylab, ' "*', names(y.factor)))
	ylim <- c(0, 1.5e-4)*y.factor
	
	plot.var.neutral.gene("default", xaxt="n", ylab=ylab, ylim=ylim, y.factor=y.factor)
	
	bottleneck.plot(Ndyn.all[["default"]], y=0.95*ylim[2], lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.95*ylim[2], cex=1.5)
	generation.axis()


	plot.var.neutral.gene("nobot", xaxt="n", ylab=ylab, ylim=ylim, y.factor=y.factor)
	
	selectionchange.plot(meansim.all[["nobot"]], y=0.95*ylim[2], cex=1.5)
	generation.axis()


	plot.var.neutral.gene("noselc", xaxt="n", ylab=ylab, ylim=ylim, y.factor=y.factor)
	
	bottleneck.plot(Ndyn.all[["noselc"]], y=0.95*ylim[2], lwd=2)
	generation.axis()

dev.off()
