#!/usr/bin/env Rscript

source("./common-fig.R")
source("../src/analysis_networks.R")

scenarios <- c("default","nomut","strongsel","strongbot")


ylim.inout.gainloss<-c(-25,25)
ylim.norm<-c(0,1.2)
ylim.nclust<-c(6,12)
ylim.Gdiff<-c(0,0.65)

ylab.inout.gainloss<-"nb connection"
ylab.norm<-"|reaction norm|"
ylab.nclust<-"nb clusters"
ylab.Gdiff<-"Change in G matrix"

pdf("figS10.pdf", width=2*length(scenarios), height = 2*4)
layout(matrix(1:(4*length(scenarios)), ncol=length(scenarios), byrow=FALSE))
par(mar=c(0.5, 0.5, 0.5, 0.5), oma=c(5, 4, 3, 0))

for (mysim in scenarios) {

  firstcol <- mysim == scenarios[1]

  plot.norm(mysim, show.quantiles=TRUE, ylim=ylim.norm,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.norm else "", xpd=if(firstcol) NA else FALSE, lty=1)
  bottleneck.plot(Ndyn.all[[mysim]], y=0, lwd=2)
  selectionchange.plot(meansim.all[[mysim]], y=0, cex=1.5)
	title(paste0(LETTERS[which(mysim == scenarios)], ": ",legname(mysim)), xpd=NA, line=2)
  
  plot.inout.gainloss(mysim, show.quantiles=TRUE, deltaG=deltaG, ylim=ylim.inout.gainloss,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.inout.gainloss else "", xpd=if(firstcol) NA else FALSE, lty=1)
  bottleneck.plot(Ndyn.all[[mysim]], y=24, lwd=2)
  selectionchange.plot(meansim.all[[mysim]], y=24, cex=1.5)
  
  plot.network.feature(mysim, what="clusters", ylim=ylim.nclust, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.nclust else "", xpd=if(firstcol) NA else FALSE, lty=1, col="black")
  bottleneck.plot(Ndyn.all[[mysim]], y=6, lwd=2)
  selectionchange.plot(meansim.all[[mysim]], y=6, cex=1.5)
  
  plot.Gdiff(mysim, deltaG=deltaG, ylim=ylim.Gdiff,xaxt="n", yaxt=if(firstcol) "s" else "n", ylab=if(firstcol) ylab.Gdiff else "", xpd=NA, lty=1, col="black")
  bottleneck.plot(Ndyn.all[[mysim]], y=0, lwd=2)
  selectionchange.plot(meansim.all[[mysim]], y=0, cex=1.5)
  generation.axis(mysim=mysim)
}

dev.off()
