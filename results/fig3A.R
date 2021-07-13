#!/usr/bin/env Rscript

# Figure: Molecular variance of neutral sites

source("./common-fig.R")

scenarios <- c("default", "nobot", "noselc")
expr.thresh <- 0.2

pdf("fig3A.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	ylab <- "Molecular variance"
	y.factor <- c('(""%*% 10^{-4})' = 10000)
	if (y.factor != 1) ylab <- parse(text=paste0('"', ylab, ' "*', names(y.factor)))
	ylim <- c(0, 1.1e-4)*y.factor
	
	plot.var.neutral(scenarios, xaxt="n", ylab=ylab, ylim=ylim, y.factor=y.factor, expr.thresh=expr.thresh, algorithm=neutral.algo)
	
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
	generation.axis()

	legend("topright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), bty="n", cex=cex.legend)

	subpanel("A")
dev.off()
