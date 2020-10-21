#!/usr/bin/env Rscript

# Figure G: number of modules

source("./commonfig.R")

connect.threshold <- 0.1
env <- 0.5
directed <- FALSE
use.cache <- TRUE

source("../src/analysis_networks.R")

comm.files <- function(out.files, cache.dir="../cache") {
	 mclapply(out.files, function(ff) {
		comm.cache.dir <- file.path(cache.dir, "comm")
		if (!dir.exists(comm.cache.dir))
			dir.create(comm.cache.dir)
		cf <- sub(".*cache/", "", dirname(ff))
		cf <- sub("/", "_", cf)
		comm.cache.file <-file.path(comm.cache.dir, paste0(cf, ".rds"))
		if (use.cache && file.exists(comm.cache.file)) {
			return(readRDS(comm.cache.file))
		} else {
			tt <- read.table(ff, header=TRUE)
			cc <- communities.dyn(tt, epsilon=connect.threshold, env=env, directed=directed, mc.cores=max(1, round(mc.cores/4)))
			saveRDS(cc, comm.cache.file, version=2)
			return(cc)
		}
	}, mc.cores=min(4, mc.cores))
}

out.files.default  <- list.files(pattern="out.*", path=list.dirs(out.dir.default,  full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.nobottle <- list.files(pattern="out.*", path=list.dirs(out.dir.nobottle, full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.noselc   <- list.files(pattern="out.*", path=list.dirs(out.dir.noselc,   full.names=TRUE, recursive=FALSE), full.names=TRUE)

comm.default  <- comm.files(out.files.default)
comm.nobottle <- comm.files(out.files.nobottle)
comm.noselc   <- comm.files(out.files.noselc)

col <- c(
#~ 	edge="violet", 
	walktrap="tomato", 
	fastgreedy="darkolivegreen",
	labelprop="lightblue")

mod.default  <- sapply(names(col), function(algo) colMeans(do.call(rbind, mclapply(comm.default, function(cc) sapply(cc, function(ccc) igraph::modularity(ccc[[algo]])), mc.cores=mc.cores))), simplify=FALSE)
mod.nobottle <- sapply(names(col), function(algo) colMeans(do.call(rbind, mclapply(comm.nobottle, function(cc) sapply(cc, function(ccc) igraph::modularity(ccc[[algo]])), mc.cores=mc.cores))), simplify=FALSE)
mod.noselc   <- sapply(names(col), function(algo) colMeans(do.call(rbind, mclapply(comm.noselc, function(cc) sapply(cc, function(ccc) igraph::modularity(ccc[[algo]])), mc.cores=mc.cores))), simplify=FALSE)

modn.default  <- sapply(names(col), function(algo) colMeans(do.call(rbind, mclapply(comm.default, function(cc) sapply(cc, function(ccc) length(igraph::communities(ccc[[algo]]))), mc.cores=mc.cores))), simplify=FALSE)
modn.nobottle <- sapply(names(col), function(algo) colMeans(do.call(rbind, mclapply(comm.nobottle, function(cc) sapply(cc, function(ccc) length(igraph::communities(ccc[[algo]]))), mc.cores=mc.cores))), simplify=FALSE)
modn.noselc   <- sapply(names(col), function(algo) colMeans(do.call(rbind, mclapply(comm.noselc, function(cc) sapply(cc, function(ccc) length(igraph::communities(ccc[[algo]]))), mc.cores=mc.cores))), simplify=FALSE)

leg <- c(default="Bottleneck + sel change", nobottle="Sel change (no bottleneck)", noselc="Bottleneck (no sel change)")

pdf("figG.pdf", width=5, height=10)
	
	layout(1:2)
	
	gen <- mean.sim.default[,"Gen"]
	plot(NULL, xlim=c(first.gen, max(gen)), ylim=c(0, max(unlist(modn.default), unlist(modn.nobottle), unlist(modn.noselc))), xlab="Generation", ylab="Nb modules")
	for (cc in names(col)) {
		my.modn.default  <- mov.avg(modn.default[[cc]],  as.numeric(names((modn.default[[cc]]))),  size=window.avg, min.gen=first.gen)
		my.modn.nobottle <- mov.avg(modn.nobottle[[cc]], as.numeric(names((modn.nobottle[[cc]]))), size=window.avg, min.gen=first.gen)
		my.modn.noselc   <- mov.avg(modn.noselc[[cc]],   as.numeric(names((modn.noselc[[cc]]))),   size=window.avg, min.gen=first.gen)
		lines(as.numeric(names(my.modn.default)),  my.modn.default,  col=col[cc], lty=1)
		lines(as.numeric(names(my.modn.nobottle)), my.modn.nobottle, col=col[cc], lty=2)
		lines(as.numeric(names(my.modn.noselc)),   my.modn.noselc,   col=col[cc], lty=3)
	}
	
	bottleneck.plot(Ndyn.default, y=0, lwd=2)
	selectionchange.plot(mean.sim.default, y=0, cex=1.5)
	legend("topleft", lty=1:3, col="black", legend=leg, bty="n", cex=0.8)
	legend("topright", pch=17, col=col, legend=names(col), bty="n", cex=0.8)


	plot(NULL, xlim=c(first.gen, max(gen)), ylim=c(0, max(unlist(mod.default), unlist(mod.nobottle), unlist(mod.noselc))), xlab="Generation", ylab="Modularity")
	for (cc in names(col)) {
		my.mod.default  <- mov.avg(mod.default[[cc]],  as.numeric(names((mod.default[[cc]]))),  size=window.avg, min.gen=first.gen)
		my.mod.nobottle <- mov.avg(mod.nobottle[[cc]], as.numeric(names((mod.nobottle[[cc]]))), size=window.avg, min.gen=first.gen)
		my.mod.noselc   <- mov.avg(mod.noselc[[cc]],   as.numeric(names((mod.noselc[[cc]]))),   size=window.avg, min.gen=first.gen)
		lines(as.numeric(names(my.mod.default)),  my.mod.default,  col=col[cc], lty=1)
		lines(as.numeric(names(my.mod.nobottle)), my.mod.nobottle, col=col[cc], lty=2)
		lines(as.numeric(names(my.mod.noselc)),   my.mod.noselc,   col=col[cc], lty=3)
	}
	
	bottleneck.plot(Ndyn.default, y=0, lwd=2)
	selectionchange.plot(mean.sim.default, y=0, cex=1.5)

#~ 	legend("topleft", lty=1:3, col="black", legend=leg, bty="n")
#~ 	legend("topright", pch=17, col=col, legend=names(col))

dev.off()
