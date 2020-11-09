#!/usr/bin/env Rscript

# Figure D: evolution of the reaction norm

source("commonfig.R")
source("../src/analysis_tools.R")

mean.norm.default  <- mean.norm.cache(out.dir.default, FUN.to.apply=abs, sliding=TRUE, window.size=window.norm)
mean.norm.nobottle <- mean.norm.cache(out.dir.nobottle, FUN.to.apply=abs, sliding=TRUE, window.size=window.norm)

pdf("figD.pdf", width=5, height=5)

plot(NULL, xlim=range(mean.sim.default[,"Gen"]), ylim=c(0, 1.2), xlab="Generation", ylab="|reaction norm|", xaxt="n")

yy.default.pp <- rowMeans(mean.norm.default[,selpattern.default=="pp"])
yy.default.ps <- rowMeans(mean.norm.default[,selpattern.default=="ps"])
yy.default.pn <- rowMeans(mean.norm.default[,selpattern.default=="pn"])

lines(as.numeric(rownames(mean.norm.default)), yy.default.pp, col=col.sel["p"], lty=lty.sce["default"])
lines(as.numeric(rownames(mean.norm.default)), yy.default.ps, col=col.sel["s"], lty=lty.sce["default"])
lines(as.numeric(rownames(mean.norm.default)), yy.default.pn, col=col.sel["n"], lty=lty.sce["default"])

yy.nobottle.pp <- rowMeans(mean.norm.nobottle[,selpattern.nobottle=="pp"])
yy.nobottle.ps <- rowMeans(mean.norm.nobottle[,selpattern.nobottle=="ps"])
yy.nobottle.pn <- rowMeans(mean.norm.nobottle[,selpattern.nobottle=="pn"])

lines(as.numeric(rownames(mean.norm.nobottle)), yy.nobottle.pp, col=col.sel["p"], lty=lty.sce["nobottle"])
lines(as.numeric(rownames(mean.norm.nobottle)), yy.nobottle.ps, col=col.sel["s"], lty=lty.sce["nobottle"])
lines(as.numeric(rownames(mean.norm.nobottle)), yy.nobottle.pn, col=col.sel["n"], lty=lty.sce["nobottle"])

bottleneck.plot(Ndyn.default, y=0, lwd=2)
selectionchange.plot(mean.sim.default, y=0, cex=1.5)

generation.axis()

legend("topright", pch=22, col=col.sel[c("p","n","s")], pt.bg=col.sel[c("p","n","s")], legend=c("Plastic->Plastic", "Plastic->Neutral", "Plastic->Stable"), bty="n", cex=0.8)
legend("topleft", lty=lty.sce[c("default","nobottle")], col="darkgray", legend=legname(c("default", "nobottle")), bty="n", cex=0.8)

dev.off()
