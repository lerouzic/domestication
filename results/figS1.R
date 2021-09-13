#!/usr/bin/env Rscript

source("./common-par.R")

s <- c(blue=10, red=50)
theta <- 0.5

pdf("figS1.pdf", width=panel.width, height=panel.height)
	par(mar=mar.notitle)
	
	plot(NULL, xlab="Gene expression", ylab="Fitness", xlim=c(0,1), ylim=c(0,1))
	
	for (ns in names(s)) {
		curve(exp(-s[ns]*(x-theta)^2), add=TRUE, xlim=c(0,1),lty=1, col=ns) 
	}
	
	legend("topright", lty=1, col=names(s), legend=paste0("s=", s), cex=cex.legend, bty="n")
dev.off()
