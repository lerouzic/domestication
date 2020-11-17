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
legend("bottomright", lty=1, col=c("black","blue"), legend=c("N", expression(N[e])), bty="n")

subpanel("A")


plot.fitness(mysim, xaxt="n", xlab="")
generation.axis()
bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)
subpanel("B")

dev.off()
