#!/usr/bin/env Rscript

# Dynamics of the effective pop size

source("./common-fig.R")

mysim <- "default"
col.opt <- c("white","blue")

plot.dynopt <- function(out.table, xlab="Generation", ...) {
	zz <- as.matrix(out.table[,grepl(colnames(out.table), pattern="FitOpt")])
	zz <- zz[,ncol(zz):1]
	image(x=out.table$Gen, y=1:ncol(zz), z=zz, col=colorRampPalette(col.opt)(100), zlim=c(0,1), xlim=c(first.gen, max(out.table$Gen)), xlab=xlab, yaxt="n", ylab="", ...)

	segreg <- selectionregime.detect(out.table)
	sel.before.dom <- sapply(strsplit(segreg, ""), "[",1)
	sel.after.dom  <- sapply(strsplit(segreg, ""), "[",2)
	sel.before.dom[1] <- sel.after.dom[1] <- "e" # The algorithm cannot know that the first guy is environment
	sel.after.dom <- ifelse(sel.before.dom != sel.after.dom, toupper(sel.after.dom), sel.after.dom)
	
	axis(2, at=ncol(zz):1, label=sel.before.dom, tick=FALSE, cex.axis=0.7, line=-0.5, las=1)
	axis(4, at=ncol(zz):1, label=sel.after.dom, tick=FALSE, cex.axis=0.7, line=-0.5, las=1)	
}

myfile <- list.files(pattern="out.*", path=onerep(outdir.all[[mysim]]),  full.names=TRUE, recursive=FALSE)
mytab <- read.table(myfile, header=TRUE)

pdf("figQ.pdf", width=panel.width, height=panel.height)
	par(mar=c(4, 1, 0.5, 1))
	plot.dynopt(mytab, xaxt="n" )

	generation.axis()
	bottleneck.plot(Ndyn.all[[mysim]], y=2, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=2, cex=1.5)
dev.off()
