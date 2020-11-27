#!/usr/bin/env Rscript

# Figure S3: evolution of gene expression variation during domestication

# Four panels: Default, no bottleneck, no change in selection, no selection. 

source("./common-fig.R")

pdf("figS3.pdf", width=1.5*panel.width, height=1.5*panel.height)
	layout(rbind(1:2,3:4))
	
	par(mar=c(0.1, 0.1, 0.5, 0.5), oma=c(5, 5, 0, 1))
	
	y.factor <- c('(""%*% 10^{-3})' = 1000)
	ylab <- "Expression variance"
	if (y.factor != 1) ylab <- parse(text=paste0('"',ylab, ' "*', names(y.factor)))
	ylim <- c(0, 1e-3)*y.factor
	
	sel.pat <- substr(selpattern.all[["default"]], 1, 1)
	sel.pat[sel.pat == "c"] <- "s" # No need to distinguish constant and stable? 
	
	### Panel A: default ###################################################
	
	plot.var.gene("default", what="expression", ylim=ylim, ylab=ylab, y.factor=y.factor, xlab="", xaxt="n", xpd=NA)
	bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
	selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
	legend("topright", lty=c(0, 1, 1, 1), col=c(0, col.sel[unique(sel.pat)]), legend=c("Before dom:", "Stable","Plastic", "Non-selected"), cex=cex.legend, bty="n")
	
	subpanel("A")
	
	### Panel B: no bottleneck #############################################
	
	plot.var.gene("nobot", what="expression", ylim=ylim, ylab="", y.factor=y.factor, xlab="", yaxt="n", xaxt="n", xpd=NA)
	selectionchange.plot(meansim.all[["nobot"]], y=0, cex=1.5)
	legend("topright", lty=c(0, lty.sel[unique(sel.pat)]), col=c(0, 1, 1, 1), legend=c("After dom:", "Stable","Plastic", "Non-selected"), cex=cex.legend, bty="n")
	
	subpanel("B")
	
	### Panel C: no selection change #######################################
	
	plot.var.gene("noselc", what="expression", ylim=ylim, ylab=ylab, y.factor=y.factor, xlab="Generation", xaxt="n", xpd=NA)
	generation.axis()
	bottleneck.plot(Ndyn.all[["noselc"]], y=0, lwd=2)
	
	subpanel("C")
	
	### Panel D: no selection ##############################################
	
	plot.var.gene("nosel", what="expression", ylim=ylim, ylab="", y.factor=y.factor, xlab="Generation", xaxt="n", yaxt="n", xpd=NA)
	generation.axis()
	bottleneck.plot(Ndyn.all[["nosel"]], y=0, lwd=2)
	subpanel("D")

dev.off()
