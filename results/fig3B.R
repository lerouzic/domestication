#!/usr/bin/env Rscript

source("../src/analysis_tools.R")

source("./common-fig.R")

scenarios <- c("default", "nobot", "noselc")

pdf(file="fig3B.pdf", height=panel.height, width=panel.width)
	par(mar=mar.notitle)
	
	y.factor <- c('(""%*% 10^{-3})' = 1000)
	ylab <- "Expression variance"
	if (y.factor != 1) ylab <- parse(text=paste0('"', ylab, ' "*', names(y.factor)))
	ylim <- c(0, 0.7e-3)*y.factor
	
	plot.var.pheno(scenarios, what="expression", y.factor=y.factor, ylab=ylab, ylim=ylim, xaxt="n", max.reps=10)
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
	
	legend(x="topright",legend = legname(scenarios), lty=lty.sce[scenarios], col=col.sce[scenarios], bty="n", cex=cex.legend)
	subpanel("B")

dev.off()
