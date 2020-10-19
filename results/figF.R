#!/usr/bin/env Rscript

# Figure F: in and out connections before/after domestication

library(parallel)
mc.cores <- min(12, detectCores()-1)

source("../src/analysis_tools.R")
source("../src/analysis_networks.R")

connect.threshold <- 0.1
env <- 0.5

my.mean.sim <- function(x) mean.sim(x, max.reps=if (detectCores() > 23) Inf else 5, mc.cores=mc.cores) # For tests 
inoutgen <- function(files, gen) 
	mclapply(files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		W <- tt[tt[,"Gen"] == gen, grepl(colnames(tt), pattern="MeanAll")]
		rm(tt); gc()
		W <- matrix(unlist(W), ncol=sqrt(length(W)), byrow=TRUE)
		ans <- try(inout.connections(W, env=env, epsilon=connect.threshold))
		if (class(ans) == "try-error") browser()
		ans
	}, mc.cores=mc.cores)
mean.inout <- function(x)
	list(connect.in = rowMeans(sapply(x, "[[", "connect.in")), connect.out=rowMeans(sapply(x, "[[", "connect.out")))

out.dir.default  <- "../cache/simDefault"
out.dir.noselc   <- "../cache/simNoselc"

out.files.default <- list.files(pattern="out.*", path=list.dirs(out.dir.default, full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.noselc  <- list.files(pattern="out.*", path=list.dirs(out.dir.noselc, full.names=TRUE, recursive=FALSE), full.names=TRUE)
if (detectCores() < 23) {
	out.files.default <- out.files.default[1:min(5, length(out.files.default))] # makes it easier for non-servers
	out.files.noselc  <- out.files.noselc[1:min(5, length(out.files.noselc))]
}
mean.sim.default  <- my.mean.sim(out.dir.default)
selpattern.default  <- selectionregime.detect(mean.sim.default)[-1]
genselchange <- selectionchange.detect(mean.sim.default)

before.dom <- inoutgen(out.files.default, gen=genselchange)
mean.before.dom <- mean.inout (before.dom)
after.dom  <- inoutgen(out.files.default, gen=mean.sim.default[nrow(mean.sim.default),"Gen"])
mean.after.dom <- mean.inout (after.dom)

# No selection change as a reference
mean.sim.noselc  <- my.mean.sim(out.dir.noselc)
selpattern.noselc <- selectionregime.detect(mean.sim.noselc)[-1]
end.sim.noselc <- inoutgen(out.files.noselc, gen=mean.sim.noselc[nrow(mean.sim.noselc),"Gen"])
mean.end.sim.noselc <- mean.inout(end.sim.noselc)


col <- c(s="blue", p="red", n="black")
leg <- c(s="Stable", p="Plastic", n="Neutral")

pdf("figF.pdf", width=5, height=5) 

xlim <- c(0, max(c(mean.before.dom$connect.in, mean.after.dom$connect.in)))
ylim <- c(0, max(c(mean.before.dom$connect.out, mean.after.dom$connect.out)))

plot(NULL, xlim=xlim, ylim=ylim, xlab="Connections to the gene", ylab="Connections from the gene")

	first.sel <- substr(selpattern.default, 1, 1)
	second.sel <- substr(selpattern.default, 2, 2)
	both.sel <- unique(selpattern.default)[-1] # removing the cc case
	points(
		x=sapply(names(col), function(nn) mean(mean.before.dom$connect.in[-1][first.sel==nn])),
		y=sapply(names(col), function(nn) mean(mean.before.dom$connect.out[-1][first.sel==nn])), 
		pch=19, col=col, cex=2)
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
		x=sapply(names(col), function(nn) mean(mean.end.sim.noselc$connect.in[-1][selpattern.noselc==nn])),
		y=sapply(names(col), function(nn) mean(mean.end.sim.noselc$connect.out[-1][selpattern.noselc==nn])), 
		pch=1, col=col, cex=2)
	
	
	legend("topleft", pch=c(15,15,15,19,17,1), col=c(col, rep("darkgray",3)), legend=c(leg, "Before domestication", "After domestication", "No domestication"))
	
dev.off()
