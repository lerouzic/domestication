#!/usr/bin/env Rscript

source("./common-fig.R")
source("../src/analysis_networks.R")

scenarios <- c("default","nomut","smallsel","strongsel","strongbot","largenet")


ylim.inout.gainloss<-c(0,20)
ylim.norm<-c(0,1.2)
ylim.nconn<-c(0,60)
ylim.Gdiff<-c(0,0.8)

ylab.inout.gainloss<-"nb connection"
ylab.norm<-"|reaction norm|"
ylab.nconn<-"nb connection"
ylab.Gdiff<-"Change in G matrix"

pdf("figS7.pdf", width=2*length(scenarios), height = 2*4)
layout(matrix(1:(4*length(scenarios)), ncol=length(scenarios), byrow=FALSE))
par(mar=c(0.5, 0.5, 0.5, 0.5), oma=c(5, 4, 3, 0))

for (mysim in scenarios) {
  firstcol <- mysim == scenarios[1]
  
  plot.norm(mysim,ylim=ylim.norm,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.norm else "", xpd=if(firstcol) NA else FALSE, lty=1)
  bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
  selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
  title(legname(mysim), xpd=NA, line=2)
  
  plot.inout.gainloss(mysim, deltaG=deltaG, ylim=ylim.inout.gainloss,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.inout.gainloss else "", xpd=if(firstcol) NA else FALSE, lty=1)
  bottleneck.plot(Ndyn.all[["default"]], y=19, lwd=2)
  selectionchange.plot(meansim.all[["default"]], y=19, cex=1.5)
  
  plot.nconn(mysim, ylim=ylim.nconn, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.nconn else "", xpd=if(firstcol) NA else FALSE, lty=1, col="black")
  bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
  selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
  
  plot.Gdiff(mysim, deltaG=deltaG, ylim=ylim.Gdiff,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.Gdiff else "", xpd=if(firstcol) NA else FALSE, lty=1, col="black")
  bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
  selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
  generation.axis()
}

dev.off()
