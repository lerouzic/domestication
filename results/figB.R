#!/usr/bin/env Rscript

# Figure B: evolution of molecular variation during domestication

# Four panels: Default, no bottleneck, no change in selection, no selection. 

source("./commonfig.R")

molec.var.default  <- molec.variation(mean.sim.default)[,-1] 
molec.var.nobottle <- molec.variation(mean.sim.nobottle)[,-1]
molec.var.noselc   <- molec.variation(mean.sim.noselc)[,-1]
molec.var.nosel    <- molec.variation(mean.sim.nosel)[,-1]


pdf("figB.pdf", width=5, height=5)
layout(rbind(1:2,3:4))

par(mar=c(0.1, 0.1, 0.1, 0.1), oma=c(5, 5, 0, 1), xpd=NA)

y.factor <- c('(""%*% 10^{-4})' = 10000)
ylab <- "Molecular variance"
if (y.factor != 1) ylab <- parse(text=paste0('"',ylab, ' "*', names(y.factor)))
ylim <- c(0, 2.5*mean(molec.var.default))*y.factor


### Panel A: default ###################################################
sel.change.gen <- selectionchange.detect(mean.sim.default)
sel.before <- substr(selpattern.default, 1, 1)

xx <- as.numeric(mean.sim.default[,"Gen"])
plot(NULL, xlim=c(first.gen, max(xx)), ylim=ylim, xlab="", ylab=ylab, xaxt="n")
for (cc in unique(sel.before)) {
	yy <- (rowMeans(molec.var.default[,sel.before==cc])*y.factor)[xx <= sel.change.gen]
	my.yy <- mov.avg(yy, xx[xx <= sel.change.gen], size=window.avg, min.gen=first.gen)
	lines(as.numeric(names(my.yy)), my.yy, lty=1, col=col.sel[cc])
}
for (cc in unique(selpattern.default)) {
	yy <- (rowMeans(molec.var.default[,selpattern.default==cc])*y.factor)[xx > sel.change.gen]
	my.yy <- mov.avg(yy, xx[xx > sel.change.gen], size=window.avg, min.gen=first.gen)
	lines(as.numeric(names(my.yy)), my.yy, lty=lty.sel[substr(cc,1,1)], col=col.sel[substr(cc,2,2)])
}
bottleneck.plot(Ndyn.default, y=0, lwd=2)
selectionchange.plot(mean.sim.default, y=0, cex=1.5)
subpanel("A")

### Panel B: no bottleneck #############################################
sel.change.gen <- selectionchange.detect(mean.sim.nobottle) # This should be the same as default
sel.before <- substr(selpattern.nobottle, 1, 1)             # idem

xx <- as.numeric(mean.sim.nobottle[,"Gen"])
plot(NULL, xlim=c(first.gen, max(xx)), ylim=ylim, xlab="", ylab="", xaxt="n", yaxt="n")
for (cc in unique(sel.before)) {
	yy <- (rowMeans(molec.var.nobottle[,sel.before==cc])*y.factor)[xx <= sel.change.gen]
	my.yy <- mov.avg(yy, xx[xx <= sel.change.gen], size=window.avg, min.gen=first.gen)
	lines(as.numeric(names(my.yy)), my.yy, lty=1, col=col.sel[cc])
}
for (cc in unique(selpattern.nobottle)) {
	yy <- (rowMeans(molec.var.nobottle[,selpattern.nobottle==cc])*y.factor)[xx > sel.change.gen]
	my.yy <- mov.avg(yy, xx[xx > sel.change.gen], size=window.avg, min.gen=first.gen)
	lines(as.numeric(names(my.yy)), my.yy, lty=lty.sel[substr(cc,1,1)], col=col.sel[substr(cc,2,2)])
}
selectionchange.plot(mean.sim.nobottle, y=0, cex=1.5)
subpanel("B")

### Panel C: no selection change #######################################
xx <- as.numeric(mean.sim.noselc[,"Gen"])
plot(NULL, xlim=c(first.gen, max(xx)), ylim=ylim, xlab="Generations", ylab=ylab, xaxt="n")
for (cc in unique(selpattern.noselc)) {
	yy <- rowMeans(molec.var.noselc[,selpattern.noselc==cc])*y.factor
	my.yy <- mov.avg(yy, xx, size=window.avg, min.gen=first.gen)
	lines(as.numeric(names(my.yy)), my.yy, lty=1, col=col.sel[cc])
}

generation.axis()
bottleneck.plot(Ndyn.noselc, y=0, lwd=2)
subpanel("C")

### Panel D: no selection ##############################################
xx <- as.numeric(mean.sim.nosel[,"Gen"])
plot(NULL, xlim=c(first.gen, max(xx)), ylim=ylim, xlab="Generations", ylab="", yaxt="n", xaxt="n")

yy <- rowMeans(molec.var.nosel)*y.factor
my.yy <- mov.avg(yy, xx, size=window.avg, min.gen=first.gen)
lines(as.numeric(names(my.yy)), my.yy, lty=1, col=col.sel["n"])

generation.axis()
bottleneck.plot(Ndyn.nosel, y=0, lwd=2)
subpanel("D")

dev.off()
