#!/usr/bin/env Rscript

source("./common-fig.R")
source("../src/analysis_networks.R")

scenarios <- c("default","step8","step24", "stab2", "nostab")

y.factor.molec <- c('(""%*% 10^{-4})' = 10000)
y.factor.expr  <- c('(""%*% 10^{-3})' = 1000)

ylab.N       <- "Population size"
ylab.fitness <- "Average fitness"
ylab.molec   <- parse(text=paste0('"', "Molecular variance", ' "*', names(y.factor.molec)))
ylab.expr    <- parse(text=paste0('"', "Expression variance", ' "*', names(y.factor.expr)))

ylab.inout.gainloss<-"nb connection"
ylab.norm<-"|reaction norm|"
ylab.nclust<-"nb clusters"
ylab.Gdiff<-"Change in G matrix"

ylim.N       <- c(0, 22000)
ylim.fitness <- c(1e-3, 1)
ylim.molec   <- c(0, y.factor.molec*1.1e-4)
ylim.expr    <- c(0, y.factor.expr*0.6e-3)

ylim.inout.gainloss<-c(-25,25)
ylim.norm<-c(0,1.2)
ylim.nclust<-c(0,14)
ylim.Gdiff<-c(0,0.8)


pdf("figS9.pdf", width=2*length(scenarios), height = 1.5*8)

layout(
	matrix(1:(8*length(scenarios)), ncol=length(scenarios), byrow=FALSE)
	)
par(mar=c(0.5, 0.5, 0.5, 0.5), oma=c(5, 4, 5, 0))

for (mysim in scenarios) {
	firstcol <- mysim == scenarios[1]
	
	plot.N(mysim, show.quantiles=TRUE, ylim=ylim.N, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.N else "", xpd=if(firstcol) NA else FALSE)
	bottleneck.plot(Ndyn.all[[mysim]], y=22000, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=22000, cex=1.5)
	title(paste0(LETTERS[which(mysim == scenarios)], ": ", legname(mysim)), xpd=NA, line=2)
	
	plot.fitness(mysim, show.quantiles=TRUE, ylim=ylim.fitness, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.fitness else "", xpd=if(firstcol) NA else FALSE, lty=1, log="y")
	bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)	

	plot.var.neutral(mysim, show.quantiles=TRUE, algorithm=neutral.algo, y.factor=y.factor.molec, ylim=ylim.molec, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.molec else "", xpd=if(firstcol) NA else FALSE)
	bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)	

	plot.var(mysim, what="expression", show.quantiles=TRUE, y.factor=y.factor.expr, ylim=ylim.expr, xlab="", xaxt="n", yaxt=if(firstcol) "s" else "n", ylab=if(firstcol) ylab.expr else "", xpd=NA)
	bottleneck.plot(Ndyn.all[[mysim]], y=0.58, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=0.58, cex=1.5)	

	plot.norm(mysim, show.quantiles=TRUE, ylim=ylim.norm,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.norm else "", xpd=if(firstcol) NA else FALSE, lty=1)
	bottleneck.plot(Ndyn.all[[mysim]], y=0, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=0, cex=1.5)
  
	plot.inout.gainloss(mysim, show.quantiles=TRUE, deltaG=deltaG, ylim=ylim.inout.gainloss,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.inout.gainloss else "", xpd=if(firstcol) NA else FALSE, lty=1)
	bottleneck.plot(Ndyn.all[[mysim]], y=24, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=24, cex=1.5)
  
	plot.network.feature(mysim, show.quantiles=TRUE, what="clusters", ylim=ylim.nclust, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.nclust else "", xpd=if(firstcol) NA else FALSE, lty=1, col="black")
	bottleneck.plot(Ndyn.all[[mysim]], y=min(ylim.nclust), lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=min(ylim.nclust), cex=1.5)
  
	plot.Gdiff(mysim, show.quantiles=TRUE, deltaG=deltaG, ylim=ylim.Gdiff,xaxt="n", yaxt=if(firstcol) "s" else "n", ylab=if(firstcol) ylab.Gdiff else "", xpd=NA, lty=1, col="black")
	bottleneck.plot(Ndyn.all[[mysim]], y=0, lwd=2)
	selectionchange.plot(meansim.all[[mysim]], y=0, cex=1.5)
	generation.axis(mysim=mysim)
}

dev.off()

