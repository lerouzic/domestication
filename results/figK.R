#!/usr/bin/env Rscript

source("../src/analysis_tools.R")

source("./common-fig.R")

scenarios <- c("default", "nobot", "noselc")

pdf(file="figK.pdf",height = 4,width = 8)
layout(t(1:2))
par(mar=c(5, 5,2 , 0.1))

###################### molecular variance plot ####################""
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

#################### expression variance plot #######################

y.factor <- c('(""%*% 10^{-3})' = 1000)
ylab <- "Expression variance"
if (y.factor != 1) ylab <- parse(text=paste0('"', ylab, ' "*', names(y.factor)))
ylim <- c(0, 0.7e-3)*y.factor

plot.var(scenarios, what="expression", y.factor=y.factor, ylab=ylab, ylim=ylim, xaxt="n")

generation.axis()
bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)

legend(x="topright",legend = legname(scenarios), lty=lty.sce[scenarios], col=col.sce[scenarios], bty="n")
subpanel("B")


dev.off()
