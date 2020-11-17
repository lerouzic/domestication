#!/usr/bin/env Rscript

# Figure E: evolution of the network complexity (number of connections)

source("./common-fig.R")

source("../src/analysis_networks.R")

scenarios <- c("default", "nobot", "noselc")

pdf("figE.pdf", width=5, height=5)

plot.nconn(scenarios, xaxt="n")

generation.axis()
bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)

legend("bottomright", lty=lty.sce[scenarios], col=col.sce[scenarios], legend=legname(scenarios), bty="n", cex=0.8)

dev.off()
