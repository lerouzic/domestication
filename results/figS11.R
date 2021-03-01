#!/usr/bin/env Rscript

source("./common-fig.R")
source("../src/analysis_networks.R")

scenarios <- c("default","smallsel","largenet","idplast","cstplast")


ylim.inout.gainloss<-c(-20,20)
ylim.norm<-c(0,1.2)
ylim.nclust<-c(6,20)
ylim.Gdiff<-c(0,0.8)

ylab.inout.gainloss<-"nb connection"
ylab.norm<-"|reaction norm|"
ylab.nclust<-"nb clusters"
ylab.Gdiff<-"Change in G matrix"

pdf("figS11.pdf", width=2*length(scenarios), height = 2*4)
layout(matrix(1:(4*length(scenarios)), ncol=length(scenarios), byrow=FALSE))
par(mar=c(0.5, 0.5, 0.5, 0.5), oma=c(5, 4, 5, 0))

mem.mc.cores <- mc.cores

for (mysim in scenarios) {

  firstcol <- mysim == scenarios[1]
  mc.cores <- mem.mc.cores
  if (mysim == "largenet") mc.cores <- min(16, mem.mc.cores)

  plot.norm(mysim,ylim=ylim.norm,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.norm else "", xpd=if(firstcol) NA else FALSE, lty=1)
  bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
  selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
	title(paste0(LETTERS[which(mysim) == scenario], ": ",legname(mysim)), xpd=NA, line=2)
  
  plot.inout.gainloss(mysim, deltaG=deltaG, ylim=ylim.inout.gainloss,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.inout.gainloss else "", xpd=if(firstcol) NA else FALSE, lty=1)
  bottleneck.plot(Ndyn.all[["default"]], y=19, lwd=2)
  selectionchange.plot(meansim.all[["default"]], y=19, cex=1.5)
  
  plot.network.feature(mysim, what="clusters", ylim=ylim.nclust, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.nclust else "", xpd=if(firstcol) NA else FALSE, lty=1, col="black")
  bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
  selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
  
  plot.Gdiff(mysim, deltaG=deltaG, ylim=ylim.Gdiff,xaxt="n", yaxt=if(firstcol) "s" else "n", ylab=if(firstcol) ylab.Gdiff else "", xpd=if(firstcol) NA else FALSE, lty=1, col="black")
  bottleneck.plot(Ndyn.all[["default"]], y=0, lwd=2)
  selectionchange.plot(meansim.all[["default"]], y=0, cex=1.5)
  generation.axis()
}

dev.off()
