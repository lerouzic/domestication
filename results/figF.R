#!/usr/bin/env Rscript

# Figure F: in and out connections before/after domestication

source("./commonfig.R")

source("../src/analysis_networks.R")

connect.threshold <- 0.1
env <- 0.5

inoutgen <- function(files, gen) {
	ans <- mclapply(files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		if (!gen %in% tt$Gen) return(NA)
		W <- tt[tt[,"Gen"] == gen, grepl(colnames(tt), pattern="MeanAll")]
		rm(tt); gc()
		W <- matrix(unlist(W), ncol=sqrt(length(W)), byrow=TRUE)
		inout.connections(W, env=env, epsilon=connect.threshold)
	}, mc.cores=mc.cores)
	ans[!sapply(ans, function(x) length(x) < 2 || is.na(x))]
}
mean.inout <- function(x)
	list(connect.in = rowMeans(sapply(x, "[[", "connect.in")), connect.out=rowMeans(sapply(x, "[[", "connect.out")))


out.files.default <- list.files(pattern="out.*", path=list.dirs(out.dir.default, full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.noselc  <- list.files(pattern="out.*", path=list.dirs(out.dir.noselc, full.names=TRUE, recursive=FALSE), full.names=TRUE)

genselchange <- selectionchange.detect(mean.sim.default)

before.dom <- inoutgen(out.files.default, gen=genselchange)
mean.before.dom <- mean.inout (before.dom)
after.dom  <- inoutgen(out.files.default, gen=mean.sim.default[nrow(mean.sim.default),"Gen"])
mean.after.dom <- mean.inout (after.dom)

# No selection change as a reference
end.sim.noselc <- inoutgen(out.files.noselc, gen=mean.sim.noselc[nrow(mean.sim.noselc),"Gen"])
mean.end.sim.noselc <- mean.inout(end.sim.noselc)

leg <- c(s="Stable", p="Plastic", n="Neutral")

pdf("figF.pdf", width=5, height=5) 

xlim <- c(0, max(c(mean.before.dom$connect.in, mean.after.dom$connect.in)))
ylim <- c(0, max(c(mean.before.dom$connect.out, mean.after.dom$connect.out)))

plot(NULL, xlim=xlim, ylim=ylim, xlab="Connections to the gene", ylab="Connections from the gene")
	
	regimes <- c("s","p","n")
	first.sel <- substr(selpattern.default, 1, 1)
	second.sel <- substr(selpattern.default, 2, 2)
	both.sel <- unique(selpattern.default)[-1] # removing the cc case
	points(
		x=sapply(regimes, function(nn) mean(mean.before.dom$connect.in[-1][first.sel==nn])),
		y=sapply(regimes, function(nn) mean(mean.before.dom$connect.out[-1][first.sel==nn])), 
		pch=19, col=col[regimes], cex=2)
	x1 <- sapply(both.sel, function(nn) mean(mean.after.dom$connect.in[-1][selpattern.default==nn]))
	y1 <- sapply(both.sel, function(nn) mean(mean.after.dom$connect.out[-1][selpattern.default==nn]))
	points(x=x1, y=y1, pch=17, col=col[substr(both.sel,2,2)], cex=2)
	arrows(
		x0=sapply(both.sel, function(nn) mean(mean.before.dom$connect.in[-1][substr(selpattern.default,1,1)==substr(nn,1,1)])),
		x1=x1,
		y0=sapply(both.sel, function(nn) mean(mean.before.dom$connect.out[-1][substr(selpattern.default,1,1)==substr(nn,1,1)])),
		y1=y1,
		col="darkgray", length=0.1)
	points(
		x=sapply(regimes, function(nn) mean(mean.end.sim.noselc$connect.in[-1][selpattern.noselc==nn])),
		y=sapply(regimes, function(nn) mean(mean.end.sim.noselc$connect.out[-1][selpattern.noselc==nn])), 
		pch=1, col=col[regimes], cex=2)
	
	legend("topleft", pch=c(15,15,15,19,17,1), col=c(col[regimes], rep("darkgray",3)), legend=c(leg, "Before domestication", "After domestication", "No domestication"))
	
dev.off()
