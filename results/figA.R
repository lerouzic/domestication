#!/usr/bin/env Rscript

# First figure: dynamics of the effective pop size and average fitness +/- sd for a random population

library(parallel)
mc.cores <- min(12, detectCores()-1)

source("../src/analysis_tools.R")

out.dir <- "../cache/simDefault"
target.reps <- NA # If NA, pick all replicates
random.rep <- NA # If NA, chose a random rep (otherwise specify it)

if (is.na(target.reps)) {
	target.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	target.reps <- target.reps[grep(basename(target.reps), pattern="rep-\\d+")]
}

if (length(target.reps %in% list.dirs(out.dir, full.names=TRUE, recursive=FALSE)) == 0) stop("Unable to find the results in ", out.dir, ".")

if (is.na(random.rep)) {
	random.rep <- sample(target.reps, 1)
}

out.files <- list.files(pattern="out.*", path=target.reps, full.names=TRUE)

stopifnot(length(out.files)>0)

tt <- mclapply(out.files, read.table, header=TRUE, mc.cores=4) # read.tables on many cores is useless, probably limited by disk speed
out.mean <- replicate.mean(tt)
out.var  <- replicate.var(tt)

mfit <- out.mean[,"MFit"]
vfit <- out.mean[,"VFit"]
gen <- out.mean[,"Gen"]

N <- get.Ndyn(random.rep)[gen]
N <- c(N[1], N) # just to deal in a bug in the naming of the parameter files

pdf("figA.pdf", width=8, height=4)

layout(t(1:2))

plot(gen, N, type="l", xlim=range(gen), ylim=c(0, max(N)), ylab="Population size", xlab="Generations")
lines(gen, N/(1+4*vfit/mfit/mfit), col="blue")
bottleneck.plot(N, y=1, lwd=2)
selectionchange.plot(out.mean, y=1, cex=1.5)
legend("bottomright", lty=1, col=c("black","blue"), legend=c("N", expression(N[e])), bty="n")

plot(NULL, xlab="Generations", ylab="Fitness", xlim=range(gen), ylim=c(0,1))
#~ lines(gen, mfit+sqrt(vfit), lty=1, col="blue")
#~ lines(gen, mfit-sqrt(vfit), lty=1, col="blue")
lines(gen, mfit, lwd=2)
bottleneck.plot(N, y=1, lwd=2)
selectionchange.plot(out.mean, y=1, cex=1.5)

dev.off()
