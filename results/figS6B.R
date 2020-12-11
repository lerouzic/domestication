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

cex.axis <- 0.5

pdf("figS6B.pdf", width=1.2*panel.width, height=panel.height)
	rl <- 0.2 #relative width of panel 1 (color scale) 
	layout(t(1:2), width=c(rl, 1-rl))
	
	par(mar=c(1, 4, 3, 0.1))
	plot.Gmat.legend(absolute=absolute)
	subpanel("B", line=1)
	
	par(mar=c(1, 2, 1, 1))
	plot.Gmat("default", c(gen.dom, gen.end), absolute=absolute, asp=1)
	# Cosmetic adjustments that depend on the exact plot dimensions
	axis(1, at=1:ng, sel.before.dom[-1], cex.axis=cex.axis, tick=FALSE, line=-2.3)
	axis(2, at=1:ng, rev(sel.before.dom[-1]), cex.axis=cex.axis, tick=FALSE, line=-1)
	axis(3, at=1:ng, sel.after.dom[-1], cex.axis=cex.axis, tick=FALSE, line=-2)
	axis(4, at=1:ng, rev(sel.after.dom[-1]), cex.axis=cex.axis, tick=FALSE, line=-1.3)

	text(2, 2, "Before domestication", pos=4)
	text(ng-2, ng-2, "Present", pos=2)



dev.off()
