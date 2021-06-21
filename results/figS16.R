#~ source("common-fig.R")
source("../src/analysis_networks.R")

# Toy example of the mutational process

col.genes <- c("green", "blue", "orange", "purple")
names.genes <- c("E", "A", "B", "C")
reg.mag <- 0.4
col.lim <- c(-reg.mag, reg.mag)

net.plot <- function(W, xlim=c(-0.2, 1.2), ylim=c(-0.2, 1.2), col=c("red", "white", "black"), text.cex=1.6, shft=0.2) {
	opar <- par(mar=c(1,1,1,1))
	stopifnot(ncol(W) == 4, nrow(W) == 4)
	rr <- colorRamp(col)
	W[W < col.lim[1]] <- col.lim[1]
	W[W > col.lim[2]] <- col.lim[2]
	
	plot(NULL, xlim=xlim, ylim=ylim, axes=FALSE, xlab="", ylab="")
	
	coords.genes <- cbind(c(0, 1, 1, 0), c(1, 1, 0, 0))
	rownames(coords.genes) <- rownames(W)
	
	text(x=coords.genes[,1], y=coords.genes[,2], rownames(coords.genes), col=col.genes, cex=text.cex)
	
	for (i in 1:nrow(W)) {
		for (j in 1:ncol(W)) {
			if (i==j) next
			cc <- coords.genes[c(i,j),]
			ww <- W[j,i]
			dd <- 0.03 * if (i > j) 1 else -1
			if (ww != 0)
				arrows( x0=cc[1,1]+shft*diff(cc[,1]) + dd, 
						x1=cc[2,1]-shft*diff(cc[,1]) + dd, 
						y0=cc[1,2]+shft*diff(cc[,2]) + dd,
						y1=cc[2,2]-shft*diff(cc[,2]) + dd,
						lwd=3, angle=if(ww < 0) 90 else 45, length=0.15,
						col=rgb(rr((ww-col.lim[1])/diff(col.lim))/255))
		}
	}
	par(mar=opar)
}

W.plot <- function(W, xlim=c(-0.2, 1.2), ylim=c(-0.2, 1.2), col=c("red", "white", "black"), cex.text=1.6, main="", bold=numeric(0), brackets=TRUE) {
	stopifnot(ncol(W) == 4, nrow(W) == 4)
	rr <- colorRamp(col)
	plot(NULL, xlim=xlim, ylim=ylim, axes=FALSE, xlab="", ylab="", main=main)
	
	for (i in 1:nrow(W)) {
		for (j in 1:ncol(W)) {
			myc <- (j-1)*ncol(W)+i
			bb <- length(bold) > 0 &&  myc %in% bold
			text(x=(j-1)/(ncol(W)-1), y=1-(i-1)/(ncol(W)-1), as.character(W[i,j]), col=rgb(rr((W[i,j]-col.lim[1])/diff(col.lim))/255), cex=0.8*cex.text, font=if(bb) 2 else 1)
		}
	}
	
	if(brackets) {
		shft1 <- 0.1
		shft2 <- 0.2
		lines(c(shft1, -shft2, -shft2, shft1), c(-shft2, -shft2, 1+shft2, 1+shft2), lwd=3)
		lines(c(1-shft1, 1+shft2, 1+shft2, 1-shft1), c(-shft2, -shft2, 1+shft2, 1+shft2), lwd=3)
	}
	
	# No way to have different axis colors without cheating
	for (i in 1:ncol(W))
		axis(2, at=rev(seq(0, 1, length.out=nrow(W)))[i], labels=rownames(W)[i], tick=FALSE, lty=0, cex.axis=cex.text, col.axis=col.genes[i], las=2)
	for (j in 1:ncol(W))
		axis(3, at=seq(0, 1, length.out=ncol(W))[j], labels=colnames(W)[j], tick=FALSE, lty=0, cex.axis=cex.text, col.axis=col.genes[j])
}

dyn.plot <- function(W, env=0.5, main="") {
	mm <- model.M2(W=W, steps=16, env=env, full=TRUE)
	
	plot(NULL, xlim=c(0, 16), ylim=c(0, 1), xlab="Time steps", ylab="Gene expression", bty="n", main=main)
	abline(h=c(0,1), col="darkgray")
	polygon(x=c(12, 16, 16, 12), y=c(0,0,1,1), border=NA, col="lightgray")

	for (i in 1:nrow(W)) {
		lines(x=0:16, mm$full[i,], col=col.genes[i], lwd=2)
	}
	
}

fit.plot <- function(thetas, s=rep(10, length(thetas)), phen=rep(NA, length(thetas)), xlab="Fitness", ylab="Gene expression", ...) {
	inv.gauss <- function(theta, s, reverse=FALSE) {
	function(y) {
	if (reverse)
		theta + sqrt(-s*log(y))/s
	else
		theta - sqrt(-s*log(y))/s
	}}
	
	opar <- par(mar=c(par("mar")[1], 0, par("mar")[3:4]))
	
	plot(NULL, xlim=c(0,1), ylim=c(0,1), xlab=xlab, ylab=ylab, bty="n", xaxt="n", ...)
	axis(1, at=c(0,1/2,1), labels=c("0", "0.5", "1"))
	abline(h=c(0,1), col="darkgray")
	for (i in seq_along(thetas)) {
		if (is.na(thetas[i]) || is.na(s[i]) || s[i]==0)
			next

		curve(inv.gauss(thetas[i], s[i])(x), add=TRUE, col=col.genes[i], lwd=2)
		curve(inv.gauss(thetas[i], s[i], TRUE)(x), add=TRUE, col=col.genes[i], lwd=2)
		if (!is.na(phen[i]))
		abline(h=phen[i], col=col.genes[i], lty=2)
	}
	par(mar=opar)
}


thetas <- c(NA, 0.7, 0.3, NA)

W.before <- W.after <- reg.mag*rbind(
	c( 0, 0, 0, 0),
    c( 1, 0,  0, 0.25),
    c(-1, 1,  0, 0),
    c( 1, 0,  0, 0))
W.after[2,4] <- W.after[2,4]*4
rownames(W.before) <- rownames(W.after) <- colnames(W.before) <- colnames(W.after) <- names.genes
      

pdf("figS16.pdf", width=12, height=6)
layout(matrix(1:8, byrow=TRUE, ncol=4), widths=c(0.2,0.2,0.4,0.1))

par(cex=0.8)

W.plot(W.before)
net.plot(W.before)
dyn.plot(W.before, main="Original network")
fit.plot(thetas, ylab="", yaxt="n", phen=model.M2(W.before)$mean)


W.plot(W.after, bold=which(W.before != W.after))
net.plot(W.after)
dyn.plot(W.after, main="Mutant network")
fit.plot(thetas, ylab="", yaxt="n", phen=model.M2(W.after)$mean)
dev.off()
