#!/usr/bin/env Rscript

# Figure D: evolution of the reaction norm

source("commonfig.R")

mean.norm.default  <- mean.norm.cache(out.dir.default, sliding=TRUE, window.size=window.avg)
mean.norm.nobottle <- mean.norm.cache(out.dir.nobottle, sliding=TRUE, window.size=window.avg)

lty <- c(default=1, nobottle=2)

pdf("figD.pdf", width=5, height=5)

plot(NULL, xlim=range(mean.sim.default[,"Gen"]), ylim=c(0, 1.2), xlab="Generation", ylab="|reaction norm|")

yy.default.pp <- rowMeans(mean.norm.default[,selpattern.default=="pp"])
yy.default.ps <- rowMeans(mean.norm.default[,selpattern.default=="ps"])
yy.default.pn <- rowMeans(mean.norm.default[,selpattern.default=="pn"])

lines(as.numeric(rownames(mean.norm.default)), yy.default.pp, col=col["p"], lty=lty["default"])
lines(as.numeric(rownames(mean.norm.default)), yy.default.ps, col=col["s"], lty=lty["default"])
lines(as.numeric(rownames(mean.norm.default)), yy.default.pn, col=col["n"], lty=lty["default"])

yy.nobottle.pp <- rowMeans(mean.norm.nobottle[,selpattern.nobottle=="pp"])
yy.nobottle.ps <- rowMeans(mean.norm.nobottle[,selpattern.nobottle=="ps"])
yy.nobottle.pn <- rowMeans(mean.norm.nobottle[,selpattern.nobottle=="pn"])

lines(as.numeric(rownames(mean.norm.nobottle)), yy.nobottle.pp, col=col["p"], lty=lty["nobottle"])
lines(as.numeric(rownames(mean.norm.nobottle)), yy.nobottle.ps, col=col["s"], lty=lty["nobottle"])
lines(as.numeric(rownames(mean.norm.nobottle)), yy.nobottle.pn, col=col["n"], lty=lty["nobottle"])

bottleneck.plot(Ndyn.default, y=0, lwd=2)
selectionchange.plot(mean.sim.default, y=0, cex=1.5)

legend("topright", pch=22, col=col, pt.bg=col, legend=c("Plastic->Plastic", "Plastic->Neutral", "Plastic->Stable"), bty="n", cex=0.8)
legend("topleft", lty=c(1,2), col="darkgray", legend=c("Bottleneck", "No bottleneck"), bty="n", cex=0.8)

dev.off()
