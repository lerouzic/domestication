#!/usr/bin/env Rscript

source("../src/analysis_tools.R")

source("./commonfig.R")

onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

molec.var.default  <- molec.variation(mean.sim.default)[,-1] 
molec.var.nobottle <- molec.variation(mean.sim.nobottle)[,-1]
molec.var.noselc   <- molec.variation(mean.sim.noselc)[,-1]

expr.var.default  <- pheno.variation(mean.sim.default)[,-1] 
expr.var.nobottle <- pheno.variation(mean.sim.nobottle)[,-1]
expr.var.noselc   <- pheno.variation(mean.sim.noselc)[,-1]

sel.change.gen <- selectionchange.detect(mean.sim.default)
selpattern<- selectionregime.detect(mean.sim.default)[-1]
sel.before <- substr(selpattern, 1, 1)
xx <- as.numeric(mean.sim.default[,"Gen"])

scenarios <- c("default", "nobottle", "noselc")

pdf(file="figK.pdf",height = 4,width = 8)
layout(t(1:2))
par(mar=c(5, 5,2 , 0.1))

###################### non selected molecular variance plot ####################""
y.factor <- c('(""%*% 10^{-4})' = 10000)
ylab <- "Molecular variance"
if (y.factor != 1){ ylab <- parse(text=paste0('"',ylab, ' "*', names(y.factor)))}
ylim <- c(0, 2.5*mean(molec.var.default))*y.factor

plot(NULL, xlim=c(first.gen, max(xx)), ylim=ylim, xlab="Generation from now", ylab=ylab ,xaxt="n")
lines(xx,rowMeans(molec.var.default [,selpattern=="nn"])*y.factor, col=col.sce["default"],  lty=lty.sce["default"])
lines(xx,rowMeans(molec.var.nobottle[,selpattern=="nn"])*y.factor, col=col.sce["nobottle"], lty=lty.sce["nobottle"])
lines(xx,rowMeans(molec.var.noselc  [,selpattern=="nn"])*y.factor, col=col.sce["noselc"],   lty=lty.sce["noselc"])

generation.axis()
bottleneck.plot(Ndyn.default, y=0, lwd=2)
selectionchange.plot(mean.sim.default, y=0, cex=1.5)

legend(x="topright",legend = legname(scenarios), lty=lty.sce[scenarios], col=col.sce[scenarios], bty="n")
subpanel("A")

#################### expression variance plot #######################

y.factor <- c('(""%*% 10^{-3})' = 1000)
ylab <- "Expression variance"
if (y.factor != 1) ylab <- parse(text=paste0('"',ylab, ' "*', names(y.factor)))
ylim <- c(0, 5*mean(expr.var.default))*y.factor

plot(NULL, xlim=c(first.gen, max(xx)), ylim=ylim, xlab="Generation from now", ylab=ylab, xaxt="n")

lines(xx,rowMeans(expr.var.default)*y.factor,  col=col.sce["default"],   lty=lty.sce["default"])
lines(xx,rowMeans(expr.var.nobottle)*y.factor, col=col.sce["nobottle"],  lty=lty.sce["nobottle"])
lines(xx,rowMeans(expr.var.noselc)*y.factor,   col=col.sce["noselc"],    lty=lty.sce["noselc"])

generation.axis()
bottleneck.plot(Ndyn.default, y=0, lwd=2)

selectionchange.plot(mean.sim.default, y=0, cex=1.5)
subpanel("B")


dev.off()
