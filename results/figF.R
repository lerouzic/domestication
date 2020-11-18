#!/usr/bin/env Rscript

# Figure F: in and out connections before/after domestication

source("./common-fig.R")

leg <- c(s="Stable", p="Plastic", n="Neutral")

pdf("figF.pdf", width=5, height=5) 
	
	regimes <- c("s","p","n")

	plot.inout.change("default", "noselc")
	
	# all pch are hard-coded... Not very clean, but the figure is used only once.
	legend("topleft", pch=c(15,15,15,19,17,1), col=c(col.sel[regimes], rep("darkgray",3)), 	
		legend=c(leg, "Before domestication", "Now", "No domestication"))
	
dev.off()
