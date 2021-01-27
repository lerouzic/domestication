#!/usr/bin/env Rscript

# Figure: Evolution of average regulation

source("./common-fig.R")

source("../src/analysis_networks.R")

scenarios <- c("default", "nobot", "noselc")

pdf("figSS1a.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.WFUN(scenarios, WFUN="mean", ylim=c(0,0.02), xaxt="n", ylab="Mean regulation effect")
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0.05, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.05, cex=1.5)
dev.off()

pdf("figSS1b.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.WFUN(scenarios, WFUN="function(w) mean(abs(w))", ylim=c(0,0.04), xaxt="n", ylab="Mean abs(regulation effect)")
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0.05, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.05, cex=1.5)
dev.off()

pdf("figSS1c.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.WFUN(scenarios, WFUN="function(w) sum( abs(w) > 0.1)", ylim=c(0,100), xaxt="n", ylab="Number of non-null connextions")
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0.05, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.05, cex=1.5)
dev.off()

pdf("figSS1d.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.WFUN(scenarios, WFUN="function(w) sum(w>0.1)", ylim=c(0,80), xaxt="n", ylab="Positive connections")
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0.05, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.05, cex=1.5)
dev.off()

pdf("figSS1e.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot.WFUN(scenarios, WFUN="function(w) sum(w < -0.1)", ylim=c(0,40), xaxt="n", ylab="Negative connections")
	
	legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), cex=cex.legend, bty="n")
	
	generation.axis()
	bottleneck.plot(Ndyn.all[["default"]], y=0.05, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0.05, cex=1.5)
dev.off()
