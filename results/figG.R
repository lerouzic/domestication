#!/usr/bin/env Rscript

# Figure G: number of modules

library(parallel)
mc.cores <- min(12, detectCores()-1)

connect.threshold <- 0.1
env <- 0.5
directed <- FALSE

source("../src/analysis_tools.R")
source("../src/analysis_networks.R")

onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]
my.mean.sim <- function(x) mean.sim(x, max.reps=if (detectCores() > 23) Inf else 5, mc.cores=mc.cores) # For tests 

comm.files <- function(out.files) {
	 mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		cc <- communities.dyn(tt, epsilon=connect.threshold, env=env, directed=directed, mc.cores=max(1, round(mc.cores/4)))
	}, mc.cores=min(4, mc.cores))
}

out.dir.default  <- "../cache/simDefault"
out.dir.nobottle <- "../cache/simNobot"
out.dir.noselc   <- "../cache/simNoselc"


out.files.default  <- list.files(pattern="out.*", path=list.dirs(out.dir.default,  full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.nobottle <- list.files(pattern="out.*", path=list.dirs(out.dir.nobottle, full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.noselc   <- list.files(pattern="out.*", path=list.dirs(out.dir.noselc,   full.names=TRUE, recursive=FALSE), full.names=TRUE)
if (detectCores() < 23) {
	out.files.default  <- out.files.default [1:min(5, length(out.files.default))] # makes it easier for non-servers
	out.files.nobottle <- out.files.nobottle[1:min(5, length(out.files.nobottle))]
	out.files.noselc   <- out.files.noselc  [1:min(5, length(out.files.noselc))]
}

comm.default  <- comm.files(out.files.default)
comm.nobottle <- comm.files(out.files.nobottle)
comm.noselc   <- comm.files(out.files.noselc)

col <- c(
#~ 	edge="violet", 
	walktrap="tomato", 
	fastgreedy="darkolivegreen",
	labelprop="lightblue")

mod.default  <- sapply(names(col), function(algo) colMeans(do.call(rbind, lapply(comm.default, function(cc) sapply(cc, function(ccc) igraph::modularity(ccc[[algo]]))))), simplify=FALSE)
mod.nobottle <- sapply(names(col), function(algo) colMeans(do.call(rbind, lapply(comm.nobottle, function(cc) sapply(cc, function(ccc) igraph::modularity(ccc[[algo]]))))), simplify=FALSE)
mod.noselc   <- sapply(names(col), function(algo) colMeans(do.call(rbind, lapply(comm.noselc, function(cc) sapply(cc, function(ccc) igraph::modularity(ccc[[algo]]))))), simplify=FALSE)

modn.default  <- sapply(names(col), function(algo) colMeans(do.call(rbind, lapply(comm.default, function(cc) sapply(cc, function(ccc) length(igraph::communities(ccc[[algo]])))))), simplify=FALSE)
modn.nobottle <- sapply(names(col), function(algo) colMeans(do.call(rbind, lapply(comm.nobottle, function(cc) sapply(cc, function(ccc) length(igraph::communities(ccc[[algo]])))))), simplify=FALSE)
modn.noselc   <- sapply(names(col), function(algo) colMeans(do.call(rbind, lapply(comm.noselc, function(cc) sapply(cc, function(ccc) length(igraph::communities(ccc[[algo]])))))), simplify=FALSE)

Ndyn.default <- get.Ndyn(onerep(out.dir.default))
mean.sim.default  <- my.mean.sim(out.dir.default)

leg <- c(default="Bottleneck + sel change", nobottle="Sel change (no bottleneck)", noselc="Bottleneck (no sel change)")


pdf("figG.pdf", width=5, height=10)
	
	layout(1:2)
	
	plot(NULL, xlim=range(mean.sim.default[,"Gen"]), ylim=c(0, max(unlist(modn.default), unlist(modn.nobottle), unlist(modn.noselc))), xlab="Generation", ylab="Nb modules")
	for (cc in names(col)) {
		lines(as.numeric(names(modn.default[[cc]])), modn.default[[cc]], col=col[cc], lty=1)
		lines(as.numeric(names(modn.nobottle[[cc]])), modn.nobottle[[cc]], col=col[cc], lty=2)
		lines(as.numeric(names(modn.noselc[[cc]])), modn.noselc[[cc]], col=col[cc], lty=3)		
	}
	
	bottleneck.plot(Ndyn.default, y=0, lwd=2)
	selectionchange.plot(mean.sim.default, y=0, cex=1.5)
	legend("topleft", lty=1:3, col="black", legend=leg, bty="n", cex=0.8)
	legend("topright", pch=17, col=col, legend=names(col), bty="n", cex=0.8)


	plot(NULL, xlim=range(mean.sim.default[,"Gen"]), ylim=c(0, max(unlist(mod.default), unlist(mod.nobottle), unlist(mod.noselc))), xlab="Generation", ylab="Modularity")
	for (cc in names(col)) {
		lines(as.numeric(names(mod.default[[cc]])), mod.default[[cc]], col=col[cc], lty=1)
		lines(as.numeric(names(mod.nobottle[[cc]])), mod.nobottle[[cc]], col=col[cc], lty=2)
		lines(as.numeric(names(mod.noselc[[cc]])), mod.noselc[[cc]], col=col[cc], lty=3)		
	}
	
	bottleneck.plot(Ndyn.default, y=0, lwd=2)
	selectionchange.plot(mean.sim.default, y=0, cex=1.5)

#~ 	legend("topleft", lty=1:3, col="black", legend=leg, bty="n")
#~ 	legend("topright", pch=17, col=col, legend=names(col))

dev.off()
