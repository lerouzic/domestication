#!/usr/bin/env Rscript

# Figure M: G matrix before/after domestication

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

cex.axis <- 0.6

pdf("figM.pdf", width=9, height=4)

layout(t(1:3), width=c(0.1, 0.45, 0.45))

par(mar=c(1, 4, 4, 0.1))

plot.Gmat.legend(absolute=absolute)

par(mar=c(0.1, 2, 4, 0.1))

plot.Gmat("default", gen.dom, absolute=absolute)
title("Before domestication")
axis(2, at=1:ng, rev(sel.before.dom[-1]), cex.axis=cex.axis, tick=FALSE, line=-1)
axis(3, at=1:ng, sel.before.dom[-1], cex.axis=cex.axis, tick=FALSE, line=-1)


plot.Gmat("default", gen.end, absolute=absolute)
title("Now")
axis(2, at=1:ng, rev(sel.after.dom[-1]), cex.axis=cex.axis, tick=FALSE, line=-1)
axis(3, at=1:ng, sel.after.dom[-1], cex.axis=cex.axis, tick=FALSE, line=-1)

dev.off()
