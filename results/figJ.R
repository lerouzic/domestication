#!/usr/bin/env Rscript

# Figure J: evolution of the number of gained/lost connections

source("./common-fig.R")

source("../src/analysis_networks.R")

scenarios <- c("default", "nobot", "noselc")

pdf("figJ.pdf", width=5, height=5)

plot.inout.gainloss(scenarios, deltaG=deltaG)

legend("topright", lty=c(1,1,lty.sce[scenarios]), col=c(col.gl, rep("black", length(scenarios))), legend=legname(c(names(col.gl), scenarios)))

bottleneck.plot(Ndyn.all[["default"]], y=30, lwd=2)
selectionchange.plot(meansim.all[["default"]], y=30, cex=1.5)

dev.off()
