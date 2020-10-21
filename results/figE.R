#!/usr/bin/env Rscript

# Figure E: evolution of the network complexity (number of connections)

# Computing the number of connections for each time point is time consuming. 
# However, a cache system (in ../cache/distW) makes further calls faster.

source("./commonfig.R")

source("../src/analysis_networks.R")

connect.threshold <- 0.1

my.mean.connect <- function(x ) {
	cache.file <- paste0(file.path(cache.dir, basename(x)), "-connect.rds")
	if (use.cache && file.exists(cache.file)) {
		return(readRDS(cache.file))
	} else {
		ans <- mean.connect(x, env=0.5, epsilon=connect.threshold, max.reps=Inf, mc.cores=mc.cores)
		saveRDS(ans, cache.file, version=2)
		return(ans)
	}
}

mean.connect.default  <- my.mean.connect(out.dir.default)
mean.connect.nobottle <- my.mean.connect(out.dir.nobottle)
mean.connect.noselc   <- my.mean.connect(out.dir.noselc)

my.mean.connect.default <- mov.avg(mean.connect.default, as.numeric(names(mean.connect.default)), size=window.avg, min.gen=first.gen)
my.mean.connect.nobottle <- mov.avg(mean.connect.nobottle, as.numeric(names(mean.connect.nobottle)), size=window.avg, min.gen=first.gen)
my.mean.connect.noselc <- mov.avg(mean.connect.noselc, as.numeric(names(mean.connect.noselc)), size=window.avg, min.gen=first.gen)

col <- c(default="black", nobottle="darkgreen", noselc="orange")
leg <- c(default="Bottleneck + sel change", nobottle="Sel change (no bottleneck)", noselc="Bottleneck (no sel change)")

pdf("figE.pdf", width=5, height=5)

plot(NULL, xlim=c(first.gen, max(mean.sim.default[,"Gen"])), ylim=c(0, max(c(mean.connect.default, mean.connect.nobottle, mean.connect.noselc))), xlab="Generation", ylab="Nb connections")

lines(as.numeric(names(my.mean.connect.default)),  my.mean.connect.default,  col=col["default"])
lines(as.numeric(names(my.mean.connect.nobottle)), my.mean.connect.nobottle, col=col["nobottle"])
lines(as.numeric(names(my.mean.connect.noselc)),   my.mean.connect.noselc,   col=col["noselc"])

bottleneck.plot(Ndyn.default, y=0, lwd=2)
selectionchange.plot(mean.sim.default, y=0, cex=1.5)

legend("topleft", pch=22, col=col, pt.bg=col, legend=leg, bty="n", cex=0.8)

dev.off()
