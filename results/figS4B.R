#!/usr/bin/env Rscript

# Figure S4B: in and out connections before/after domestication

source("./common-fig.R")

leg <- c(s="Stable", p="Plastic", n="Neutral")

pdf("figS4B.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	regimes <- c("s","p","n")

	plot.inout.change("default", "noselc")
	
	# all pch are hard-coded... Not very clean, but the figure is used only once.
	legend("topleft", pch=c(15,15,15,19,17,1), col=c(col.sel[regimes], rep("darkgray",3)), 	
		legend=c(leg, "Before domestication", "Now", "No domestication"), cex=cex.legend, bty="n")
		
	subpanel("B")
	
dev.off()
