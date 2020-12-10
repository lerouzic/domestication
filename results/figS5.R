#!/usr/bin/env Rscript

# Figure S5: G matrix before/after domestication

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

pdf("figS5.pdf", width=2.2*panel.width, height=panel.height)
	rl <- 0.45 #relative width of panels 2 and 3 
	layout(t(1:3), width=c(0.08, rl, rl))
	
	par(mar=c(1, 4, 1, 0.1))
	plot.Gmat.legend(absolute=absolute)
	subpanel("A", adj=0.5, col="white")
	
	par(mar=c(1, 2, 1, 1))
	plot.Gmat("default", c(gen.dom, gen.end), absolute=absolute, asp=1)
	axis(1, at=1:ng, sel.before.dom[-1], cex.axis=cex.axis, tick=FALSE, line=-1)
	axis(2, at=1:ng, rev(sel.before.dom[-1]), cex.axis=cex.axis, tick=FALSE, line=-1)
	axis(3, at=1:ng, sel.after.dom[-1], cex.axis=cex.axis, tick=FALSE, line=-1)
	axis(4, at=1:ng, rev(sel.after.dom[-1]), cex.axis=cex.axis, tick=FALSE, line=-1)

	text(2, 2, "Before domestication", pos=4)
	text(ng-2, ng-2, "Present", pos=2)

	my.col.sel <- col.sel[unique(tolower(c(sel.before.dom,sel.after.dom)))]
	my.col.sel <- lighten.col(my.col.sel[!is.na(my.col.sel)], factor=0.3)

	par(mar=mar.notitle)
	
	plot.Gmat.all("default", gens=c(gen.dom, gen.end), absolute=absolute, 
		cols=my.col.sel, sel1=sel.before.dom, sel2=tolower(sel.after.dom), ylim=c(0,0.8))
	subpanel("B")
		
	par(fig=c(1-rl*0.2,1-rl*0.05,0.8,0.95), new=TRUE, mar=c(0,0,0,0))
	plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab="", ylab="", axes=FALSE, asp=1)
	rasterImage(as.raster(outer(my.col.sel, my.col.sel, FUN=avgcol)),xleft=0,xright=1,ybottom=0,ytop=1,interpolate=FALSE)
	axis(3, at=(0.5+0:(length(my.col.sel)-1))/length(my.col.sel), names(my.col.sel), tick=FALSE, line=-1)
	axis(2, at=(0.5+0:(length(my.col.sel)-1))/length(my.col.sel), rev(names(my.col.sel)), tick=FALSE, line=-1)

dev.off()
