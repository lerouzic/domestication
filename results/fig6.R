#!/usr/bin/env Rscript

source("./common-fig.R")
source("../src/analysis_networks.R")

scenarios <- c("default","africe","pemil","tomato")

y.factor.molec <- c('(""%*% 10^{-4})' = 10000)
y.factor.expr  <- c('(""%*% 10^{-3})' = 1000)

ylab.N       <- "Population size"
ylab.molec   <- parse(text=paste0('"', "Mol var ", ' "*', names(y.factor.molec)))
ylab.expr    <- parse(text=paste0('"', "Expr var", ' "*', names(y.factor.expr)))
ylab.norm    <- "|reaction norm|"
ylab.inout.gainloss <- "Nb connections"
ylab.Gcor   <- "Genetic correlation"

ylim.N       <- c(0, 22000)
ylim.molec   <- c(0, y.factor.molec*1.1e-4)
ylim.expr    <- c(0, y.factor.expr*0.6e-3)
ylim.norm    <- c(0,1.2)
ylim.inout.gainloss <- c(-25,25)
ylim.Gcor   <- c(0,0.6)

pdf("fig6.pdf", width=2*length(scenarios), height = 2*4)

layout(matrix(1:(6*length(scenarios)), ncol=length(scenarios), byrow=FALSE))
par(mar=c(0.5, 0.5, 0.5, 0.5), oma=c(4, 4, 3, 0))

for (mysim in scenarios) {
	firstcol <- mysim == scenarios[1]
	
	plot.N(mysim, show.quantiles=TRUE, ylim=ylim.N, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.N else "", xpd=if(firstcol) NA else FALSE)
	title(paste0(LETTERS[which(mysim == scenarios)], ": ", if (mysim=="default") "Maize" else legname(mysim)), xpd=NA, line=2)
	
	plot.var.neutral(mysim, show.quantiles=TRUE, algorithm=neutral.algo, y.factor=y.factor.molec, ylim=ylim.molec, xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.molec else "", xpd=if(firstcol) NA else FALSE)
	bottleneck.plot(Ndyn.all[[mysim]], lwd=1, style="vlines", col="darkgray")
	
	plot.var(mysim, what="expression", show.quantiles=TRUE, y.factor=y.factor.expr, ylim=ylim.expr, xaxt="n", xlab="", yaxt=if(firstcol) "s" else "n", ylab=if(firstcol) ylab.expr else "", xpd=NA)
	bottleneck.plot(Ndyn.all[[mysim]], lwd=1, style="vlines", col="darkgray")
	
	plot.norm(mysim, show.quantiles=TRUE, ylim=ylim.norm,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.norm else "", xpd=if(firstcol) NA else FALSE, lty=1)
	bottleneck.plot(Ndyn.all[[mysim]], lwd=1, style="vlines", col="darkgray")
	
	plot.inout.gainloss(mysim, show.quantiles=TRUE, deltaG=deltaG, ylim=ylim.inout.gainloss,xaxt="n", yaxt=if(firstcol) "s" else "n", xlab="", ylab=if(firstcol) ylab.inout.gainloss else "", xpd=if(firstcol) NA else FALSE, lty=1)
	bottleneck.plot(Ndyn.all[[mysim]], lwd=1, style="vlines", col="darkgray")
	
	plot.Gcor(mysim, show.quantiles=TRUE, ylim=ylim.Gcor,xaxt="n", yaxt=if(firstcol) "s" else "n", ylab=if(firstcol) ylab.Gcor else "", xpd=NA, lty=1, col="black")
	bottleneck.plot(Ndyn.all[[mysim]], lwd=1, style="vlines", col="darkgray")
	generation.axis(mysim=mysim)
}

dev.off()
