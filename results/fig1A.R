#!/usr/bin/env Rscript

source("../src/analysis_tools.R")

source("./common-fig.R")

scenarios <- c("default", "nobot", "noselc")

pdf(file="fig1A.pdf", height=panel.height, width=panel.width)
layout(t(1:2))
par(mar=c(5, 5,2 , 0.1))

y.factor <- c('(""%*% 10^{-4})' = 10000)

ylab <- "Molecular variance"
y.factor <- c('(""%*% 10^{-4})' = 10000)
if (y.factor != 1) ylab <- parse(text=paste0('"', ylab, ' "*', names(y.factor)))
ylim <- c(0, 1e-4)*y.factor

plot.var(scenarios, what="molecular", y.factor=y.factor, ylab=ylab, ylim=ylim, xaxt="n")

generation.axis()
bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)

subpanel("A")

dev.off()
