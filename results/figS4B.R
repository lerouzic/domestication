#!/usr/bin/env Rscript

# Figure S4B: in and out connections before/after domestication

source("./common-fig.R")

leg <- c(s="Stable", p="Plastic", n="Neutral")

pdf("figS4B.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	regimes <- c("s","p","n")

	plot.inout.change("default", "noselc", xlim=c(0, 6), ylim=c(0,6))
	
	# all pch are hard-coded... Not very clean, but the figure is used only once.
	legend("topright", pch=c(15,15,15,19,17,1), col=c(col.sel[regimes], rep("darkgray",3)), 	
		legend=c(leg, "Before domestication", "Now", legname("noselc")), cex=cex.legend, bty="n")
		
	subpanel("B")
	
dev.off()
