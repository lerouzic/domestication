#!/usr/bin/env Rscript

# First figure: dynamics of the effective pop size and average fitness +/- sd for a random population

source("./commonfig.R")

mfit <- mean.sim.default[,"MFit"]
vfit <- mean.sim.default[,"VFit"]
gen <-  mean.sim.default[,"Gen"]

N <- Ndyn.default
names(N) <- as.character(as.numeric(names(N)))

my.N     <- mov.avg(N[gen+1], gen, size=window.avg, min.gen=first.gen)
my.mfit  <- mov.avg(mfit, gen, size=window.avg, min.gen=first.gen)
my.vfit  <- mov.avg(vfit, gen, size=window.avg, min.gen=first.gen)

pdf("figA.pdf", width=8, height=4)

layout(t(1:2))

plot(as.numeric(names(my.N)), my.N, type="l", xlim=range(as.numeric(names(my.N))), ylim=c(0, max(N)), ylab="Population size", xlab="Generations") # no moving average
lines(as.numeric(names(my.N)), my.N/(1+4*my.vfit/my.mfit/my.mfit), col="blue")
bottleneck.plot(my.N, y=1, lwd=2)
selectionchange.plot(mean.sim.default, y=1, cex=1.5)
legend("bottomright", lty=1, col=c("black","blue"), legend=c("N", expression(N[e])), bty="n")

plot(NULL, xlab="Generations", ylab="Fitness", xlim=range(gen), ylim=c(0,1))
#~ lines(gen, mfit+sqrt(vfit), lty=1, col="blue")
#~ lines(gen, mfit-sqrt(vfit), lty=1, col="blue")
lines(as.numeric(names(my.mfit)), my.mfit, lwd=2)
bottleneck.plot(N, y=1, lwd=2)
selectionchange.plot(mean.sim.default, y=1, cex=1.5)

dev.off()
