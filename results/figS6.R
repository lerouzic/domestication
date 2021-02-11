#!/usr/bin/env Rscript

# Figure H: gene networks before/after domestication

source("./common-fig.R")

source("../src/analysis_networks.R")

mysim <- "default"

genselchange <- selectionchange.detect(meansim.all[[mysim]])

Wlist.before.dom <- Wgen.files.cache(outdir.all[[mysim]], gen=genselchange)
Wlist.after.dom  <- Wgen.files.cache(outdir.all[[mysim]], gen=meansim.all[[mysim]][nrow(meansim.all[[mysim]]),"Gen"])

segreg <- selectionregime.detect(meansim.all[[mysim]])
sel.before.dom <- sapply(strsplit(segreg, ""), "[",1)
sel.after.dom <- sapply(strsplit(segreg, ""), "[",2)
sel.before.dom[sel.before.dom == "c"] <- "s" # constant and stable should be the same
sel.after.dom[sel.after.dom == "c"] <- "s"
sel.before.dom[1] <- sel.after.dom[1] <- "e" # The algorithm cannot know that the first guy is environment

numconn.before.dom <- mean.numconn.groups(Wlist.before.dom, sel.before.dom, mc.cores=mc.cores)
numconn.after.dom  <- mean.numconn.groups(Wlist.after.dom,  sel.after.dom,  mc.cores=mc.cores)

# Complex name tags (including the number of genes)
cn <- colnames(numconn.before.dom$plus)
group.names.before.dom <- setNames(cn, paste0(toupper(cn), "[", table(sel.before.dom)[cn], "]"))
group.names.after.dom  <- setNames(cn, paste0(toupper(cn), "[", table(sel.after.dom)[cn],  "]"))

pdf("figS6.pdf", width=2*panel.width, height=panel.height)
	layout(t(1:2))
	
	plot.numconn.groups(numconn.before.dom, group.names=group.names.before.dom, thresh=0.05)
	title("Before domestication")
	
	plot.numconn.groups(numconn.after.dom, numconn.ref=numconn.before.dom, group.names=group.names.after.dom, thresh=0.02)
	title("Present")
	
dev.off()
