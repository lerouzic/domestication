#!/usr/bin/env Rscript

# Figure B: evolution of molecular variation during domestication

# Four panels: Default, no bottleneck, no change in selection, no selection. 

library(parallel)
mc.cores <- min(12, detectCores()-1)

source("../src/analysis_tools.R")

mean.sim <- function(out.dir, max.reps=5) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)

	tt <- mclapply(out.files[1:(if (detectCores() > 24 || length(out.files) < max.reps) length(out.files) else max.reps)], read.table, header=TRUE, mc.cores=4) # read.tables on many cores is useless, probably limited by disk speed
	ans <- replicate.mean(tt)
	rm(tt)
	gc()
	return(ans)
}

onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

out.dir.default  <- "../cache/simDefault"
out.dir.nobottle <- "../cache/simNobot"
out.dir.noselc   <- "../cache/simNoselc"
out.dir.nosel    <- "../cache/simNosel"

mean.sim.default  <- mean.sim(out.dir.default)
mean.sim.nobottle <- mean.sim(out.dir.nobottle)
mean.sim.noselc   <- mean.sim(out.dir.noselc)
mean.sim.nosel    <- mean.sim(out.dir.nosel)

selpattern.default  <- selectionregime.detect(mean.sim.default)[-1] # The first gene is the environmental signal
selpattern.nobottle <- selectionregime.detect(mean.sim.nobottle)[-1]
selpattern.noselc   <- selectionregime.detect(mean.sim.noselc)[-1]
selpattern.nosel    <- selectionregime.detect(mean.sim.nosel)[-1]

molec.var.default  <- molec.variation(mean.sim.default)[,-1] 
molec.var.nobottle <- molec.variation(mean.sim.nobottle)[,-1]
molec.var.noselc   <- molec.variation(mean.sim.noselc)[,-1]
molec.var.nosel    <- molec.variation(mean.sim.nosel)[,-1]

Ndyn.default    <- get.Ndyn(onerep(out.dir.default))
Ndyn.nobottle   <- get.Ndyn(onerep(out.dir.nobottle))
Ndyn.noselc     <- get.Ndyn(onerep(out.dir.noselc))
Ndyn.nosel      <- get.Ndyn(onerep(out.dir.nosel))

col <- c(c="black", s="gray", p="red", n="darkgreen")
lty <- c(c=1, s=1, p=2, n=3)

pdf("figB.pdf", width=5, height=5)
layout(rbind(1:2,3:4))

par(mar=c(0.1, 0.1, 0.1, 0.1), oma=c(5, 5, 0, 0), xpd=NA)

ylim <- c(0, 4*mean(molec.var.default))

### Panel A: default ###################################################
sel.change.gen <- selectionchange.detect(mean.sim.default)
sel.before <- substr(selpattern.default, 1, 1)

xx <- as.numeric(mean.sim.default[,"Gen"])
plot(NULL, xlim=range(xx), ylim=ylim, xlab="", ylab="Molecular variance", xaxt="n")
for (cc in unique(sel.before)) {
	yy <- rowMeans(molec.var.default[,sel.before==cc])
	lines(xx[xx <= sel.change.gen], yy[xx <= sel.change.gen], lty=1, col=col[cc])
}
for (cc in unique(selpattern.default)) {
	yy <- rowMeans(molec.var.default[,selpattern.default==cc])
	lines(xx[xx > sel.change.gen], yy[xx > sel.change.gen], lty=lty[substr(cc,1,1)], col=col[substr(cc,2,2)])
}
bottleneck.plot(Ndyn.default, y=0, lwd=2)
selectionchange.plot(mean.sim.default, y=0, cex=1.5)

### Panel B: no bottleneck #############################################
sel.change.gen <- selectionchange.detect(mean.sim.nobottle) # This should be the same as default
sel.before <- substr(selpattern.nobottle, 1, 1)             # idem

xx <- as.numeric(mean.sim.nobottle[,"Gen"])
plot(NULL, xlim=range(xx), ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n")
for (cc in unique(sel.before)) {
	yy <- rowMeans(molec.var.nobottle[,sel.before==cc])
	lines(xx[xx <= sel.change.gen], yy[xx <= sel.change.gen], lty=1, col=col[cc])
}
for (cc in unique(selpattern.nobottle)) {
	yy <- rowMeans(molec.var.nobottle[,selpattern.nobottle==cc])
	lines(xx[xx > sel.change.gen], yy[xx > sel.change.gen], lty=lty[substr(cc,1,1)], col=col[substr(cc,2,2)])
}
selectionchange.plot(mean.sim.nobottle, y=0, cex=1.5)

### Panel C: no selection change #######################################
xx <- as.numeric(mean.sim.noselc[,"Gen"])
plot(NULL, xlim=range(xx), ylim=ylim, xlab="Generations", ylab="Molecular variance")
for (cc in unique(selpattern.noselc)) {
	yy <- rowMeans(molec.var.noselc[,selpattern.noselc==cc])
	lines(xx, yy, lty=1, col=col[cc])
}
bottleneck.plot(Ndyn.noselc, y=0, lwd=2)

### Panel D: no selection ##############################################
xx <- as.numeric(mean.sim.nosel[,"Gen"])
plot(NULL, xlim=range(xx), ylim=ylim, xlab="Generations", ylab="", yaxt="n")
axis(1)

yy <- rowMeans(molec.var.noselc)
lines(xx, yy, lty=1, col=col["c"])

bottleneck.plot(Ndyn.nosel, y=0, lwd=2)

dev.off()
