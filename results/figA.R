#!/usr/bin/env Rscript

# First figure: dynamics of the effective pop size and average fitness +/- sd for a random population

source("./commonfig.R")

mfit <- mean.sim.default[,"MFit"]
vfit <- mean.sim.default[,"VFit"]
gen <-  mean.sim.default[,"Gen"]

pdf("figA.pdf", width=8, height=4)

layout(t(1:2))

plot(as.numeric(names(Ndyn.default)), Ndyn.default, type="l", xlim=c(first.gen, max(gen)), ylim=c(0, max(Ndyn.default)), xaxt="n", ylab="Population size", xlab="Generations") # no moving average
lines(gen, Ndyn.default[as.character(gen)]/(1+4*vfit/mfit/mfit), col="blue")

generation.axis()
bottleneck.plot(Ndyn.default, y=1, lwd=2)
selectionchange.plot(mean.sim.default, y=1, cex=1.5)

subpanel("A")
legend("bottomright", lty=1, col=c("black","blue"), legend=c("N", expression(N[e])), bty="n")

plot(NULL, xlab="Generations", ylab="Fitness", xlim=c(0, max(gen)), ylim=c(0,1), xaxt="n")
#~ lines(gen, mfit+sqrt(vfit), lty=1, col="blue")
#~ lines(gen, mfit-sqrt(vfit), lty=1, col="blue")
lines(gen, mfit)

generation.axis()
bottleneck.plot(Ndyn.default, y=1, lwd=2)
selectionchange.plot(mean.sim.default, y=1, cex=1.5)
subpanel("B")

dev.off()
