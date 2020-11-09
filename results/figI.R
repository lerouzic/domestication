#!/usr/bin/env Rscript

# Figure I: coexpression networks before/after domestication

source("./commonfig.R")

source("../src/analysis_networks.R")

cutoff <- 0.1
covtransf <- cov2cor # could be identity

out.files.default <- list.files(pattern="out.*", path=list.dirs(out.dir.default, full.names=TRUE, recursive=FALSE), full.names=TRUE)
genselchange <- selectionchange.detect(mean.sim.default)

Clist.before.dom <- lapply(Ggen.cache(out.files.default, gen=genselchange), covtransf)
Clist.after.dom  <- lapply(Ggen.cache(out.files.default, gen=mean.sim.default[nrow(mean.sim.default),"Gen"]), covtransf)

segreg <- selectionregime.detect(mean.sim.default)
sel.before.dom <- sapply(strsplit(segreg, ""), "[",1)
sel.after.dom  <- sapply(strsplit(segreg, ""), "[",2)
sel.before.dom[sel.before.dom == "c"] <- "s" # constant and stable should be the same
sel.after.dom [sel.after.dom  == "c"] <- "s"
sel.before.dom[1] <- sel.after.dom[1] <- "e" # The algorithm cannot know that the first guy is environment

numcorr.before.dom    <- mean.numcorrgen.groups(Clist.before.dom, sel.before.dom, cutoff=cutoff, mc.cores=mc.cores)
numcorr.justafter.dom <- mean.numcorrgen.groups(Clist.before.dom, sel.after.dom,  cutoff=cutoff, mc.cores=mc.cores)
numcorr.after.dom     <- mean.numcorrgen.groups(Clist.after.dom,  sel.after.dom,  cutoff=cutoff, mc.cores=mc.cores)

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
