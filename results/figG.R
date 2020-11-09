#!/usr/bin/env Rscript

# Figure G: number of modules

source("./commonfig.R")

connect.threshold <- 0.1
env <- 0.5
directed <- FALSE
use.cache <- TRUE

source("../src/analysis_networks.R")

comm.files <- function(out.files) {
	 mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		cc <- communities.dyn.cache(tt, epsilon=connect.threshold, env=env, directed=directed, mc.cores=max(1, round(mc.cores/4)))
		return(cc)
		}, mc.cores=min(4, mc.cores))
}

out.files.default  <- list.files(pattern="out.*", path=list.dirs(out.dir.default,  full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.nobottle <- list.files(pattern="out.*", path=list.dirs(out.dir.nobottle, full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.noselc   <- list.files(pattern="out.*", path=list.dirs(out.dir.noselc,   full.names=TRUE, recursive=FALSE), full.names=TRUE)

comm.default  <- comm.files(out.files.default)
comm.nobottle <- comm.files(out.files.nobottle)
comm.noselc   <- comm.files(out.files.noselc)

col.model <- c(
#~ 	edge="violet", 
	walktrap="tomato", 
	fastgreedy="darkolivegreen",
	labelprop="lightblue")

mod.default  <- sapply(names(col.model), function(algo) colMeans(do.call(rbind, mclapply(comm.default, function(cc) sapply(cc, function(ccc) igraph::modularity(ccc[[algo]])), mc.cores=mc.cores))), simplify=FALSE)
mod.nobottle <- sapply(names(col.model), function(algo) colMeans(do.call(rbind, mclapply(comm.nobottle, function(cc) sapply(cc, function(ccc) igraph::modularity(ccc[[algo]])), mc.cores=mc.cores))), simplify=FALSE)
mod.noselc   <- sapply(names(col.model), function(algo) colMeans(do.call(rbind, mclapply(comm.noselc, function(cc) sapply(cc, function(ccc) igraph::modularity(ccc[[algo]])), mc.cores=mc.cores))), simplify=FALSE)

modn.default  <- sapply(names(col.model), function(algo) colMeans(do.call(rbind, mclapply(comm.default, function(cc) sapply(cc, function(ccc) length(igraph::communities(ccc[[algo]]))), mc.cores=mc.cores))), simplify=FALSE)
modn.nobottle <- sapply(names(col.model), function(algo) colMeans(do.call(rbind, mclapply(comm.nobottle, function(cc) sapply(cc, function(ccc) length(igraph::communities(ccc[[algo]]))), mc.cores=mc.cores))), simplify=FALSE)
modn.noselc   <- sapply(names(col.model), function(algo) colMeans(do.call(rbind, mclapply(comm.noselc, function(cc) sapply(cc, function(ccc) length(igraph::communities(ccc[[algo]]))), mc.cores=mc.cores))), simplify=FALSE)

scenarios <- c("default","nobottle","noselc")

pdf("figG.pdf", width=5, height=10)
	
	layout(1:2)
	
	gen <- mean.sim.default[,"Gen"]
	plot(NULL, xlim=c(first.gen, max(gen)), ylim=c(0, max(unlist(modn.default), unlist(modn.nobottle), unlist(modn.noselc))), xlab="Generation", ylab="Nb modules", xaxt="n")
	for (cc in names(col.model)) {
		my.modn.default  <- mov.avg(modn.default[[cc]],  as.numeric(names((modn.default[[cc]]))),  size=window.avg, min.gen=first.gen)
		my.modn.nobottle <- mov.avg(modn.nobottle[[cc]], as.numeric(names((modn.nobottle[[cc]]))), size=window.avg, min.gen=first.gen)
		my.modn.noselc   <- mov.avg(modn.noselc[[cc]],   as.numeric(names((modn.noselc[[cc]]))),   size=window.avg, min.gen=first.gen)
		lines(as.numeric(names(my.modn.default)),  my.modn.default,  col=col.model[cc], lty=lty.sce["default"])
		lines(as.numeric(names(my.modn.nobottle)), my.modn.nobottle, col=col.model[cc], lty=lty.sce["nobottle"])
		lines(as.numeric(names(my.modn.noselc)),   my.modn.noselc,   col=col.model[cc], lty=lty.sce["noselc"])
	}
	
	generation.axis()
	bottleneck.plot(Ndyn.default, y=0, lwd=2)
	selectionchange.plot(mean.sim.default, y=0, cex=1.5)
	legend("topleft", lty=lty.sce[scenarios], col="black", legend=legname(scenarios), bty="n", cex=0.8)
	legend("topright", pch=17, col=col.model, legend=names(col.model), bty="n", cex=0.8)
	
	subpanel("A")


	plot(NULL, xlim=c(first.gen, max(gen)), ylim=c(0, max(unlist(mod.default), unlist(mod.nobottle), unlist(mod.noselc))), xlab="Generation", ylab="Modularity", xaxt="n")
	for (cc in names(col.model)) {
		my.mod.default  <- mov.avg(mod.default[[cc]],  as.numeric(names((mod.default[[cc]]))),  size=window.avg, min.gen=first.gen)
		my.mod.nobottle <- mov.avg(mod.nobottle[[cc]], as.numeric(names((mod.nobottle[[cc]]))), size=window.avg, min.gen=first.gen)
		my.mod.noselc   <- mov.avg(mod.noselc[[cc]],   as.numeric(names((mod.noselc[[cc]]))),   size=window.avg, min.gen=first.gen)
		lines(as.numeric(names(my.mod.default)),  my.mod.default,  col=col.model[cc], lty=lty.sce["default"])
		lines(as.numeric(names(my.mod.nobottle)), my.mod.nobottle, col=col.model[cc], lty=lty.sce["nobottle"])
		lines(as.numeric(names(my.mod.noselc)),   my.mod.noselc,   col=col.model[cc], lty=lty.sce["noselc"])
	}
	
	generation.axis()
	bottleneck.plot(Ndyn.default, y=0, lwd=2)
	selectionchange.plot(mean.sim.default, y=0, cex=1.5)

	subpanel("B")

#~ 	legend("topleft", lty=1:3, col="black", legend=leg, bty="n")
#~ 	legend("topright", pch=17, col=col, legend=names(col))

dev.off()
