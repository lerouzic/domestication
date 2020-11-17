#!/usr/bin/env Rscript

# Figure I: coexpression networks before/after domestication

source("./common-fig.R")

source("../src/analysis_networks.R")

mysim <- "default"

covtransf <- cov2cor # could be identity

out.files <- list.files(pattern="out.*", path=list.dirs(outdir.all[[mysim]], full.names=TRUE, recursive=FALSE), full.names=TRUE)
genselchange <- selectionchange.detect(meansim.all[[mysim]])

Clist.before.dom <- lapply(Ggen.files.cache(out.files, gen=genselchange), covtransf)
Clist.after.dom  <- lapply(Ggen.files.cache(out.files, gen=meansim.all[[mysim]][nrow(meansim.all[[mysim]]),"Gen"]), covtransf)

segreg <- selectionregime.detect(meansim.all[[mysim]])
sel.before.dom <- sapply(strsplit(segreg, ""), "[",1)
sel.after.dom  <- sapply(strsplit(segreg, ""), "[",2)
sel.before.dom[sel.before.dom == "c"] <- "s" # constant and stable should be the same
sel.after.dom [sel.after.dom  == "c"] <- "s"
sel.before.dom[1] <- sel.after.dom[1] <- "e" # The algorithm cannot know that the first guy is environment

numcorr.before.dom    <- mean.numcorrgen.groups(Clist.before.dom, sel.before.dom, mc.cores=mc.cores)
numcorr.after.dom     <- mean.numcorrgen.groups(Clist.after.dom,  sel.after.dom,  mc.cores=mc.cores)

cn <- colnames(numcorr.before.dom$plus)
group.names.before.dom <- setNames(cn, paste0(toupper(cn), "[", table(sel.before.dom)[cn], "]"))
group.names.after.dom  <- setNames(cn, paste0(toupper(cn), "[", table(sel.after.dom)[cn],  "]"))

pdf("figI.pdf", width=8, height=4)
	layout(t(1:2))
	
	plot.numconn.groups(numcorr.before.dom, group.names=group.names.before.dom, directed=FALSE, pos.shift.minus=0.5)
	title("Before domestication")
	plot.numconn.groups(numcorr.after.dom, numconn.ref=numcorr.before.dom, group.names=group.names.after.dom, directed=FALSE, pos.shift.minus=0.5)
	title("Now")
dev.off()
