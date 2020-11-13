#!/usr/bin/env Rscript

# Figure M: G matrix before/after domestication

source("./commonfig.R")

source("../src/analysis_networks.R")

# Plotting correlations
covtransf <- cov2cor # vs identity for covariances, but the color scale will not match

col.cor <- c("red", "white", "green") # for correlations -1, 0, 1

meanG <- function(listG, transf=identity) {
	arr <- do.call(abind, c(lapply(listG, transf), list(along=3)))
	rowMeans(arr, dims=2)
}

out.files.default <- list.files(pattern="out.*", path=list.dirs(out.dir.default, full.names=TRUE, recursive=FALSE), full.names=TRUE)

Glist.before.dom <- lapply(Ggen.cache(out.files.default, gen=selectionchange.detect(mean.sim.default)), covtransf)
Glist.after.dom  <- lapply(Ggen.cache(out.files.default, gen=mean.sim.default[nrow(mean.sim.default),"Gen"]), covtransf)

Gmean.before.dom <- meanG(Glist.before.dom)
Gmean.after.dom  <- meanG(Glist.after.dom)

segreg <- selectionregime.detect(mean.sim.default)

sel.before.dom <- sapply(strsplit(segreg, ""), "[",1)
sel.after.dom  <- sapply(strsplit(segreg, ""), "[",2)
sel.before.dom[sel.before.dom == "c"] <- "s" # constant and stable should be the same
sel.after.dom [sel.after.dom  == "c"] <- "s"
sel.before.dom[1] <- sel.after.dom[1] <- "e" # The algorithm cannot know that the first guy is environment
sel.after.dom <- ifelse(sel.before.dom != sel.after.dom, toupper(sel.after.dom), sel.after.dom)

Gmean.before.toplot <- Gmean.before.dom
Gmean.after.toplot  <- Gmean.after.dom

ng <- ncol(Gmean.before.toplot)-1
cols <- colorRampPalette(col.cor)(1001)
cex.axis <- 0.6

pdf("figM.pdf", width=9, height=4)

layout(t(1:3), width=c(0.1, 0.45, 0.45))

par(mar=c(1, 4, 4, 0.1))

image(y=seq(-1, 1, length.out=length(cols)), z=t(as.matrix(seq(-1, 1, length.out=length(cols)-1))), zlim=c(-1,1), col=cols, xaxt="n", xlab="", ylab="Genetic correlation")

par(mar=c(0.1, 2, 4, 0.1))

image(x=1:ng, y=1:ng, t(Gmean.before.toplot[(ng+1):2,2:(ng+1)]), axes=FALSE, xlab="", ylab="", zlim=c(-1, 1), col=cols)
title("Before domestication")
axis(2, at=1:ng, rev(sel.before.dom[-1]), cex.axis=cex.axis, tick=FALSE, line=-1)
axis(3, at=1:ng, sel.before.dom[-1], cex.axis=cex.axis, tick=FALSE, line=-1)


image(x=1:ng, y=1:ng, t(Gmean.after.toplot[(ng+1):2,2:(ng+1)]), axes=FALSE, xlab="", ylab="", zlim=c(-1, 1), col=cols)
title("After domestication")
axis(2, at=1:ng, rev(sel.after.dom[-1]), cex.axis=cex.axis, tick=FALSE, line=-1)
axis(3, at=1:ng, sel.after.dom[-1], cex.axis=cex.axis, tick=FALSE, line=-1)

dev.off()
