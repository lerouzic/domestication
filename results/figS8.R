#!/usr/bin/env Rscript

source("./common-fig.R")

scenarios <- c("default","nomut","strongsel","strongbot")

y.factor.molec <- c('(""%*% 10^{-4})' = 10000)
y.factor.expr  <- c('(""%*% 10^{-3})' = 1000)

ylab.N       <- "Population size"
ylab.fitness <- "Average fitness"
ylab.molec   <- parse(text=paste0('"', "Molecular variance", ' "*', names(y.factor.molec)))
ylab.expr    <- parse(text=paste0('"', "Expression variance", ' "*', names(y.factor.expr)))

ylim.N       <- c(0, 22000)
ylim.fitness <- c(1e-3, 1)
ylim.molec   <- c(0, y.factor.molec*1.1e-4)
ylim.expr    <- c(0, y.factor.expr*1e-3)

pdf("figS8.pdf", width=2*length(scenarios), height = 2*4)

layout(matrix(1:(4*length(scenarios)), ncol=length(scenarios), byrow=FALSE))
par(mar=c(0.5, 0.5, 0.5, 0.5), oma=c(5, 4, 3, 0))

for (mysim in scenarios) {
	firstcol <- mysim == scenarios[1]
	
	plot.N(mysim, ylim=ylim.N, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.N else "", xpd=if(firstcol) NA else FALSE)
	bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)
	title(legname(mysim), xpd=NA, line=2)
	
	plot.fitness(mysim, ylim=ylim.fitness, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.fitness else "", xpd=if(firstcol) NA else FALSE, lty=1, log="y")
	bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)	

	plot.var.gene(mysim, what="molecular", y.factor=y.factor.molec, ylim=ylim.molec, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.molec else "", xpd=if(firstcol) NA else FALSE)
	bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)	

	plot.var.gene(mysim, what="expression", y.factor=y.factor.expr, ylim=ylim.expr, xaxt="n", yaxt=if(firstcol) "s" else "n", ylab=if(firstcol) ylab.expr else "", xpd=NA)
	bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)	
	generation.axis()	
}

dev.off()

