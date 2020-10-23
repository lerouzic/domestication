#!/usr/bin/env Rscript

# Figure H: gene networks before/after domestication

source("./commonfig.R")

source("../src/analysis_networks.R")

connect.threshold <- 0.1
env <- 0.5

Wgen <- function(files, gen) {
	mclapply(files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		W <- tt[tt[,"Gen"] == gen, grepl(colnames(tt), pattern="MeanAll")]
		rm(tt); gc()
		W <- matrix(unlist(W), ncol=sqrt(length(W)), byrow=TRUE)
	}, mc.cores=1)
}

Wgen.cache <- function(files, gen) {
	cache.fun(Wgen, files=files, gen=gen, cache.subdir="Wgen")
}

out.files.default <- list.files(pattern="out.*", path=list.dirs(out.dir.default, full.names=TRUE, recursive=FALSE), full.names=TRUE)
genselchange <- selectionchange.detect(mean.sim.default)

Wlist.before.dom <- Wgen.cache(out.files.default, gen=genselchange)
Wlist.after.dom  <- Wgen.cache(out.files.default, gen=mean.sim.default[nrow(mean.sim.default),"Gen"])

segreg <- selectionregime.detect(mean.sim.default)
sel.before.dom <- sapply(strsplit(segreg, ""), "[",1)
sel.after.dom <- sapply(strsplit(segreg, ""), "[",2)
sel.before.dom[sel.before.dom == "c"] <- "s" # constant and stable should be the same
sel.after.dom[sel.after.dom == "c"] <- "s"
sel.before.dom[1] <- sel.after.dom[1] <- "e" # The algorithm cannot know that the first guy is environment

numconn.before.dom <- mean.numconn.groups(Wlist.before.dom, sel.before.dom)
numconn.justafter.dom <- mean.numconn.groups(Wlist.before.dom, sel.after.dom)
numconn.after.dom <- mean.numconn.groups(Wlist.after.dom, sel.after.dom)


pdf("figH.pdf", width=12, height=4)
layout(t(1:3))

plot.numconn.groups(numconn.before.dom)
title("Before domestication")
plot.numconn.groups(numconn.justafter.dom)
title("Immediately after")
plot.numconn.groups(numconn.after.dom)
title("Now")

dev.off()
