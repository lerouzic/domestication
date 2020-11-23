#!/usr/bin/env Rscript

# First figure: dynamics of the effective pop size and average fitness +/- sd for a random population

source("./common-fig.R")

mysim <- "default"

pdf("figA.pdf", width=8, height=4)

layout(t(1:2))

plot.N(mysim, xaxt="n", xlab="")
generation.axis()
bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)
legend("bottomright", lty=1, col=c("black","blue"), legend=c("N", expression(N[e])), cex=cex.legend)

subpanel("A")

scenarios <- c("default", "nobot", "noselc")

plot.fitness(scenarios, xaxt="n", xlab="")
generation.axis()
bottleneck.plot(Ndyn.all[["default"]], y=1, lwd=2)
selectionchange.plot(meansim.all[["default"]], y=1, cex=1.5)
legend("bottomright", lty=lty.sce[scenarios], col="black", legend=legname(scenarios), cex=cex.legend)
subpanel("B")

dev.off()
