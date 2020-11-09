#!/usr/bin/env Rscript

# Figure L: Evolution of genetic variance-covariance

library(ade4) # for mantel.rtest

source("./commonfig.R")

source("../src/analysis_networks.R")

delta.gen <- 1000

cor2dist <- function(r) as.dist(sqrt(2*(1-r)))
covtransf <- function(x) cor2dist(cov2cor(x))

Mantel.obs <- function(lG) {
	ans <-matrix(NA, nrow=length(lG)-1, ncol=length(lG[[1]]))
	for (i in 1:(length(lG)-1)) {
		for (j in 1:length(lG[[i]])) {
			ans[i, j] <- mantel.rtest(lG[[i]][[j]], lG[[i+1]][[j]], nrepet=1)$obs
		}
	}
	ans
}

PC.partvariance <- function(lG, whichPC=1) {
	ans <- matrix(NA, nrow=length(lG), ncol=length(lG[[1]]))
	for (i in 1:length(lG)) {
		for (j in 1:length(lG[[1]])) {
			ee <- eigen(lG[[i]][[j]])
			ans[i,j] <- ee$values[whichPC]/sum(ee$values)
		}
	}
	ans
}

out.files.default  <- list.files(pattern="out.*", path=list.dirs(out.dir.default,  full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.nobottle <- list.files(pattern="out.*", path=list.dirs(out.dir.nobottle, full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.noselc   <- list.files(pattern="out.*", path=list.dirs(out.dir.noselc,   full.names=TRUE, recursive=FALSE), full.names=TRUE)

genselchange <- selectionchange.detect(mean.sim.default)

gens <- mean.sim.default[,"Gen"]
mygens <- gens[seq(1, length(gens), length.out=round(gens[length(gens)] / delta.gen))]

allG.default   <- lapply(mygens, function(g) Ggen.cache(out.files.default,  g))
allG.nobottle  <- lapply(mygens, function(g) Ggen.cache(out.files.nobottle, g))
allG.noselc    <- lapply(mygens, function(g) Ggen.cache(out.files.noselc,   g))

mt.default     <- Mantel.obs(lapply(allG.default, function(x) lapply(x, covtransf)))
mt.nobottle    <- Mantel.obs(lapply(allG.nobottle, function(x) lapply(x, covtransf)))
mt.noselc      <- Mantel.obs(lapply(allG.noselc, function(x) lapply(x, covtransf)))

pv.default     <- PC.partvariance(lapply(allG.default, function(x) lapply(x, cov2cor)))
pv.nobottle    <- PC.partvariance(lapply(allG.nobottle, function(x) lapply(x, cov2cor)))
pv.noselc      <- PC.partvariance(lapply(allG.noselc, function(x) lapply(x, cov2cor)))

scenarios <- c("default","nobottle","noselc")

pdf("figL.pdf", width=8, height=4)
	layout(rbind(1:2))
	
	plot(NULL, xlim=c(first.gen, max(gens)), ylim=c(0,1), xlab="Generations", ylab="Change in G matrix", xaxt="n")
	lines(mygens[-1], 1-rowMeans(mt.default), col=col.sce["default"], lty=lty.sce["default"])
	lines(mygens[-1], 1-rowMeans(mt.nobottle), col=col.sce["nobottle"], lty=lty.sce["nobottle"])
	lines(mygens[-1], 1-rowMeans(mt.noselc), col=col.sce["noselc"], lty=lty.sce["noselc"])
	
	generation.axis()
	bottleneck.plot(Ndyn.default, y=0, lwd=2)
	selectionchange.plot(mean.sim.default, y=0, cex=1.5)
	
	legend(x="topright",legend = legname(scenarios), lty=lty.sce[scenarios], col=col.sce[scenarios], bty="n")
	subpanel("A")
	
	
	plot(NULL, xlim=c(first.gen, max(gens)), ylim=c(0,0.4), xlab="Generations", ylab="Proportion PC1/total variance", xaxt="n")
	lines(mygens, rowMeans(pv.default), col=col.sce["default"], lty=lty.sce["default"])
	lines(mygens, rowMeans(pv.nobottle), col=col.sce["nobottle"], lty=lty.sce["nobottle"])
	lines(mygens, rowMeans(pv.noselc), col=col.sce["noselc"], lty=lty.sce["noselc"])
	
	generation.axis()
	bottleneck.plot(Ndyn.default, y=0, lwd=2)
	selectionchange.plot(mean.sim.default, y=0, cex=1.5)

	subpanel("B")

dev.off()



