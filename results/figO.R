#!/usr/bin/env Rscript

source("./commonsupfig.R")

pdf("figO.pdf",width=4*length(simtoplot),height = 4*3)
layout(matrix(1:(3*length(simtoplot)),ncol=length(simtoplot),nrow =3,byrow=F))

for (tag in simtag){
  ##### plot type fig A ######
  mfit<-eval(parse(text=paste0("mean.sim.",tag)))[,"MFit"]
  vfit<-eval(parse(text=paste0("mean.sim.",tag)))[,"VFit"]
  gen<-eval(parse(text=paste0("mean.sim.",tag)))[,"Gen"]
  N<-eval(parse(text=paste0("Ndyn.",tag)))
  
  names(N) <- as.character(as.numeric(names(N)))
  my.N     <- mov.avg(N[gen+1], gen, size=window.avg, min.gen=first.gen)
  my.mfit  <- mov.avg(mfit, gen, size=window.avg, min.gen=first.gen)
  my.vfit  <- mov.avg(vfit, gen, size=window.avg, min.gen=first.gen)
  if(tag==simtag[1]){par(mar=c(0.1,5,4,0.1))
    plot(NULL,xlab = "", ylab="Fitness", xlim=range(gen), ylim=c(0,1),main=tag,xaxt="n")
  } else {par(mar=c(0.1,0.1,4,0.1))
    plot(NULL,xlab = "", ylab="", xlim=range(gen), ylim=c(0,1),xaxt="n",main = tag,yaxt="n")  
  }
  lines(as.numeric(names(my.mfit)), my.mfit, lwd=2)
  bottleneck.plot(N, y=1, lwd=2)
  selectionchange.plot(mean.sim.data, y=1, cex=1.5)
  
  ##### plot type fig B ########
  sel.change.gen <- selectionchange.detect(mean.sim.data)
  sel.before <- substr(selpattern.default, 1, 1)
  xx <- as.numeric(mean.sim.data[,"Gen"])
  
  y.factor <- c('(""%*% 10^{-4})' = 10000)
  ylab <- "Molecular variance"
  molec.var.default  <- molec.variation(mean.sim.default)[,-1]
  if (y.factor != 1){ ylab <- parse(text=paste0('"',ylab, ' "*', names(y.factor)))}
  ylim <- c(0, 3*mean(molec.var.default))*y.factor
  
  if(tag==simtag[1]){par(mar=c(0.1,5,0.1,0.1))
    plot(NULL, xlim=c(first.gen, max(xx)), ylim=ylim, xlab="", ylab=ylab, xaxt="n")
  } else {par(mar=c(0.1,0.1,0.1,0.1))
    plot(NULL, xlim=c(first.gen, max(xx)), ylim=ylim, xlab="", ylab="", xaxt="n",yaxt="n")
  }
  molec.var<-molec.variation(eval(parse(text=paste0("mean.sim.",tag))))
  
  for (cc in unique(sel.before)) {
    yy <- (rowMeans(molec.var[,which(sel.before==cc),drop=F])*y.factor)[xx <= sel.change.gen]
    my.yy <- mov.avg(yy, xx[xx <= sel.change.gen], size=window.avg, min.gen=first.gen)
    lines(as.numeric(names(my.yy)), my.yy, lty=1, col=col.sel[cc])
  }
  for (cc in unique(selpattern)) {
    yy <- (rowMeans(molec.var[,which(selpattern==cc),drop=F])*y.factor)[xx > sel.change.gen]
    my.yy <- mov.avg(yy, xx[xx > sel.change.gen], size=window.avg, min.gen=first.gen)
    lines(as.numeric(names(my.yy)), my.yy, lty=lty.sel[substr(cc,1,1)], col=col.sel[substr(cc,2,2)])
  }
  bottleneck.plot(Ndyn, y=0, lwd=2)
  selectionchange.plot(mean.sim.data, y=0, cex=1.5)
  
  ######## type fig C plot ###########
  expr.var.default  <- pheno.variation(mean.sim.default)[,-1]
  y.factor <- c('(""%*% 10^{-3})' = 1000)
  ylab <- "Expression variance"
  if (y.factor != 1) ylab <- parse(text=paste0('"',ylab, ' "*', names(y.factor)))
  ylim <- c(0, 5*mean(expr.var.default))*y.factor
  
  expr.var<-pheno.variation(eval(parse(text=paste0("mean.sim.",tag))))
  
  if(tag==simtag[1]){par(mar=c(0.1,5,0.1,0.1))
    plot(NULL, xlim=range(xx), ylim=ylim, xlab="", ylab=ylab,xaxt="n")
  } else {par(mar=c(0.1,0.1,0.1,0.1))
    plot(NULL, xlim=range(xx), ylim=ylim, xlab="", ylab="",yaxt="n",xaxt="n")
  }
  
  for (cc in unique(sel.before)) {
    yy <- (rowMeans(expr.var[,which(sel.before==cc),drop=F])*y.factor)[xx <= sel.change.gen]
    my.yy <- mov.avg(yy, xx[xx <= sel.change.gen], size=window.avg, min.gen=first.gen)
    lines(xx[xx <= sel.change.gen], yy[xx <= sel.change.gen], lty=1, col=col.sel[cc])
  }
  for (cc in unique(selpattern)) {
    yy <- (rowMeans(expr.var[,which(selpattern==cc),drop=F])*y.factor)[xx > sel.change.gen]
    my.yy <- mov.avg(yy, xx[xx > sel.change.gen], size=window.avg, min.gen=first.gen)
    lines(as.numeric(names(my.yy)), my.yy, lty=lty.sel[substr(cc,1,1)], col=col.sel[substr(cc,2,2)])
  }
  bottleneck.plot(Ndyn, y=0, lwd=2)
  selectionchange.plot(mean.sim.data, y=0, cex=1.5)
  
}

dev.off()

