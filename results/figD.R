#!/usr/bin/env Rscript

# Figure D: evolution of the reaction norm

source("common-fig.R")

pdf("figD.pdf", width=5, height=5)

plot.norm(c("default", "nobot"), xlab="Generation", ylab="|reaction norm|", xaxt="n")

bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
generation.axis()

legend("topright", pch=22, col=col.sel[c("p","n","s")], pt.bg=col.sel[c("p","n","s")], legend=c("Plastic->Plastic", "Plastic->Neutral", "Plastic->Stable"), bty="n", cex=0.8)
legend("topleft", lty=lty.sce[c("default","nobot")], col="darkgray", legend=legname(c("default", "nobot")), bty="n", cex=0.8)

dev.off()
