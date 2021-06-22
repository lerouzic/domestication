#!/usr/bin/env Rscript

# Dynamics of the effective pop size

source("./common-fig.R")

mysims <- c("default", "africe", "pemil", "tomato")

pdf("figS1.pdf", width=length(mysims)*panel.width, height=panel.height)
	layout(t(seq_along(mysims)))
	par(mar=c(mar.title), cex=1)
	
	for (mysim in mysims) {
		first <- mysim == mysims[1]
		plot.N(mysim, xaxt="n", xlab="Generation", ylim=c(0,22000), ylab=if(first) "Population size" else "", main=if (mysim == "default") "Maize" else legname(mysim))

		generation.axis(mysim=mysim)
		bottleneck.plot(Ndyn.all[[mysim]], y=1, lwd=2)
		selectionchange.plot(meansim.all[[mysim]], y=1, cex=1.5)
		if (first)
			legend("bottomright", lty=1, col=c("black","blue"), legend=c("N", expression(N[e])), cex=cex.legend, bty="n")
	}
	
dev.off()
