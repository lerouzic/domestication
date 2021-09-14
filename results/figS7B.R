#!/usr/bin/env Rscript

# Figure: in and out connections before/after domestication

source("./common-fig.R")

leg <- c(s="Stable", p="Plastic", n="Neutral")

pdf("figS7B.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	regimes <- c("s","p","n")

	plot.inout.change("default", xlim=c(0, 6), ylim=c(0,4))
	
	# all pch are hard-coded... Not very clean, but the figure is used only once.
	legend("bottomright", pch=c(15,15,15,19,17), col=c(col.sel[regimes], rep("darkgray",2)), 	
	legend=c(leg, "Before domestication", "Present"), cex=cex.legend, bty="n")
		
	subpanel("B")
	
dev.off()
