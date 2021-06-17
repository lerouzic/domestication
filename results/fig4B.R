#!/usr/bin/env Rscript

# Figure: G matrix before/after domestication

source("./common-fig.R")

source("../src/analysis_networks.R")

# Plotting correlations
absolute <- TRUE     # absolute value for correlations? 
mysim <- "default"

gen.dom <- selectionchange.detect(meansim.all[[mysim]])
gen.end <- meansim.all[[mysim]][nrow(meansim.all[[mysim]]),"Gen"]
segreg <- selectionregime.detect(meansim.all[[mysim]])

ng <- length(segreg)-1

sel.before.dom <- sapply(strsplit(segreg, ""), "[",1)
sel.after.dom  <- sapply(strsplit(segreg, ""), "[",2)
sel.before.dom[sel.before.dom == "c"] <- "s" # constant and stable should be the same
sel.after.dom [sel.after.dom  == "c"] <- "s"
sel.before.dom[1] <- sel.after.dom[1] <- "e" # The algorithm cannot know that the first guy is environment
sel.after.dom <- ifelse(sel.before.dom != sel.after.dom, toupper(sel.after.dom), sel.after.dom)

cex.axis <- 0.7

plot.axis <- function(side, x, ...) {
	# R does not like narrow labels for axes. This is a workaround,
	# calling 'axis' for odd and even labels successively
	odd  <- seq(1, length(x), 2)
	even <- seq(2, length(x), 2)
	axis(side=side, at=(1:ng)[odd], label=x[odd], cex.axis=cex.axis, tick=FALSE, ...)
	axis(side=side, at=(1:ng)[even], label=x[even], cex.axis=cex.axis, tick=FALSE, ...)
}

pdf("fig4B.pdf", width=1.2*panel.width, height=panel.height)
	rl <- 0.2 #relative width of panel 1 (color scale) 
	layout(t(1:2), width=c(rl, 1-rl))
	
	par(mar=c(1, 4, mar.notitle[3], 0.1))
	plot.Gmat.legend(absolute=absolute)
	subpanel("B", line=1)
	
	par(mar=c(1, 1, mar.notitle[3], 1))
	plot.Gmat("default", c(gen.dom, gen.end), absolute=absolute, asp=1)
	# Cosmetic adjustments that depend on the exact plot dimensions
	plot.axis(1, sel.before.dom[-1], line=-1.7)
	plot.axis(2, rev(sel.before.dom[-1]), line=-1)
	plot.axis(3, sel.after.dom[-1], line=-1.3)
	plot.axis(4, rev(sel.after.dom[-1]), line=-1.3)	

	text(2, 2, "Before domestication", pos=4)
	text(ng-2, ng-2, "Present", pos=2)



dev.off()
