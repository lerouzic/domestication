#!/usr/bin/env Rscript

# Figure J: evolution of the number of gained/lost connections

source("./commonfig.R")

source("../src/analysis_networks.R")

delta.gen <- 1000
col.gl <- c(gain="orange", loss="darkgray")
lty.gl <- c(default=1, nobottle=2, noselc=3)

Wgen <- function(files, gen) {
	ans <- mclapply(files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		if (!gen %in% tt[,"Gen"]) return(NA)
		W <- tt[tt[,"Gen"] == gen, grepl(colnames(tt), pattern="MeanAll")]
		rm(tt); gc()
		W <- matrix(unlist(W), ncol=sqrt(length(W)), byrow=TRUE)
	}, mc.cores=mc.cores)
	names(ans) <- files
	ans[!sapply(ans, function(x) length(x) == 1 && is.na(x))]
}

Wgen.cache <- function(files, gen) {
	cache.fun(Wgen, files=files, gen=gen, cache.subdir="Wgen")
}

gainorloss <- function(lW) {
	.gains <- function(m1, m2) sum(m1 == 0 & m2 != 0)
	.losses <- function(m1, m2) sum(m1 != 0 & m2 == 0)
	ans <- list(
		gain=matrix(NA, nrow=length(lW)-1, ncol=length(lW[[1]])),
		loss=matrix(NA, nrow=length(lW)-1, ncol=length(lW[[1]])))
	for (i in 1:(length(lW)-1)) {
		for (j in 1:length(lW[[i]])) {
			ans$gain[i, j] <- .gains(lW[[i]][[j]], lW[[i+1]][[j]])
			ans$loss[i, j] <- .losses(lW[[i]][[j]], lW[[i+1]][[j]])
		}
	}
	ans
}

out.files.default  <- list.files(pattern="out.*", path=list.dirs(out.dir.default,  full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.nobottle <- list.files(pattern="out.*", path=list.dirs(out.dir.nobottle, full.names=TRUE, recursive=FALSE), full.names=TRUE)
out.files.noselc   <- list.files(pattern="out.*", path=list.dirs(out.dir.noselc,   full.names=TRUE, recursive=FALSE), full.names=TRUE)

gens <- mean.sim.default[,"Gen"]
mygens <- gens[seq(1, length(gens), length.out=round(gens[length(gens)] / delta.gen))]

allW.default   <- setNames(lapply(mygens, function(g) Wgen.cache(out.files.default,  g)), as.character(mygens))
allW.nobottle  <- setNames(lapply(mygens, function(g) Wgen.cache(out.files.nobottle, g)), as.character(mygens))
allW.noselc    <- setNames(lapply(mygens, function(g) Wgen.cache(out.files.noselc,   g)), as.character(mygens))
allcW.default  <- lapply(allW.default, function(gw)  lapply(gw, cleanW.cache))
allcW.nobottle <- lapply(allW.nobottle, function(gw) lapply(gw, cleanW.cache))
allcW.noselc   <- lapply(allW.noselc, function(gw)   lapply(gw, cleanW.cache))
gl.default     <- gainorloss(allcW.default)
gl.nobottle    <- gainorloss(allcW.nobottle)
gl.noselc      <- gainorloss(allcW.noselc)


pdf("figJ.pdf", width=5, height=5)

plot(NULL, xlim=c(first.gen, max(mygens)), ylim=range(c(unlist(gl.default), unlist(gl.nobottle), unlist(gl.noselc))),  xlab="Generations", ylab="Number connections")

lines(mygens[-1], rowMeans(gl.default$gain),  col=col.gl["gain"], lty=lty.gl["default"])
lines(mygens[-1], rowMeans(gl.nobottle$gain), col=col.gl["gain"], lty=lty.gl["nobottle"])
lines(mygens[-1], rowMeans(gl.noselc$gain),   col=col.gl["gain"], lty=lty.gl["noselc"])

lines(mygens[-1], rowMeans(gl.default$loss),  col=col.gl["loss"], lty=lty.gl["default"])
lines(mygens[-1], rowMeans(gl.nobottle$loss), col=col.gl["loss"], lty=lty.gl["nobottle"])
lines(mygens[-1], rowMeans(gl.noselc$loss),   col=col.gl["loss"], lty=lty.gl["noselc"])

legend("topright", lty=c(1,1,lty.gl), col=c(col.gl, rep("black", length(lty.gl))), legend=legname(c(names(col.gl), names(lty.gl))))

bottleneck.plot(Ndyn.default, y=30, lwd=2)
selectionchange.plot(mean.sim.default, y=30, cex=1.5)

dev.off()
