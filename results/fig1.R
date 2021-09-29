#~ source("common-fig.R")
source("../src/analysis_networks.R")

# Toy example of the mutational process

col.genes <- c("darkgreen", "black", "black", "black")
names.genes <- c("E", "A", "B", "C")
env.default <- 0.5
env.alt     <- 0.8
col.lim <- c(-0.5, 0.5)

net.plot <- function(W, xlim=c(-0.2, 1.2), ylim=c(-0.2, 1.2), col=c("black", "white", "black"), text.cex=1.6, shft=0.2, bold=numeric(0), font.caption=rep(1,ncol(W))) {
	opar <- par(mar=c(3,1,1,1))
	stopifnot(ncol(W) == 4, nrow(W) == 4)
	rr <- colorRamp(col)
	W[W < col.lim[1]] <- col.lim[1]
	W[W > col.lim[2]] <- col.lim[2]
	
	plot(NULL, xlim=xlim, ylim=ylim, axes=FALSE, xlab="", ylab="")
	
	coords.genes <- cbind(c(0, 1, 1, 0), c(1, 1, 0, 0))
	rownames(coords.genes) <- rownames(W)
	
	text(x=coords.genes[,1], y=coords.genes[,2], rownames(coords.genes), col=col.genes, cex=text.cex, font=font.caption)
	
	for (i in 1:nrow(W)) {
		for (j in 1:ncol(W)) {
			if (i==j) next
			cc <- coords.genes[c(i,j),]
			ww <- W[j,i]
			myc <- (i-1)*ncol(W)+j
			bb <- length(bold) > 0 &&  myc %in% bold
			dd <- 0.03 * if (i > j) 1 else -1
			if (ww != 0)
				arrows( x0=cc[1,1]+shft*diff(cc[,1]) + dd, 
						x1=cc[2,1]-shft*diff(cc[,1]) + dd, 
						y0=cc[1,2]+shft*diff(cc[,2]) + dd,
						y1=cc[2,2]-shft*diff(cc[,2]) + dd,
						lwd=3, angle=if(ww < 0) 90 else 45, length=0.15,
						col=if (bb) "red" else rgb(rr((ww-col.lim[1])/diff(col.lim))/255))
		}
	}
	par(mar=opar)
}

W.plot <- function(W, xlim=c(-0.2, 1.2), ylim=c(-0.2, 1.2), col=c("black", "black", "black"), cex.text=1.6, main="", bold=numeric(0), bold.caption=numeric(0), brackets=TRUE, show.zeros=FALSE, first.line=FALSE) {
	stopifnot(ncol(W) == 4, nrow(W) == 4)
	rr <- colorRamp(col)
	plot(NULL, xlim=xlim, ylim=ylim, axes=FALSE, xlab="", ylab="", main=main)
	
	if (!first.line) {
		Wb <- matrix(FALSE, ncol=ncol(W), nrow=nrow(W))
		Wb[bold] <- TRUE
		bold <- which(Wb[-1,])
		W <- W[-1,]
	}
	
	for (i in 1:nrow(W)) {
		for (j in 1:ncol(W)) {
			myc <- (j-1)*nrow(W)+i
			bb <- length(bold) > 0 &&  myc %in% bold
			if (W[i,j] != 0 || show.zeros)
				text(
					x=(j-1)/(ncol(W)-1), 
					y=1-(i-1)/(nrow(W)-1), 
					as.character(W[i,j]), 
#~ 					col=rgb(rr((W[i,j]-col.lim[1])/diff(col.lim))/255), font=if(bb) 2 else 1
					cex=0.8*cex.text, col=if(bb) "red" else "black"
			)
		}
	}
	
	if(brackets) {
		shft1 <- 0.1
		shft2 <- 0.2
		lines(c(shft1, -shft2, -shft2, shft1), c(-shft2, -shft2, 1+shft2, 1+shft2), lwd=3)
		lines(c(1-shft1, 1+shft2, 1+shft2, 1-shft1), c(-shft2, -shft2, 1+shft2, 1+shft2), lwd=3)
	}
	
	# No way to have different axis colors without cheating
	for (i in 1:nrow(W))
		axis(2, at=rev(seq(0, 1, length.out=nrow(W)))[i], labels=rownames(W)[i], tick=FALSE, lty=0, cex.axis=cex.text, col.axis=col.genes[i+if(first.line) 0 else 1], las=2)
	for (j in 1:ncol(W))
		axis(3, at=seq(0, 1, length.out=ncol(W))[j], labels=colnames(W)[j], tick=FALSE, lty=0, cex.axis=cex.text, col.axis=col.genes[j], font.axis=if(j %in% bold.caption) 2 else 1 )
}

dyn.plot <- function(W, env=env.default, main="", gene.names=TRUE) {
	mm <- model.M2(W=W, steps=16, env=env, full=TRUE)
	
	plot(NULL, xlim=c(0, 16), ylim=c(0, 1), xlab="Time steps", ylab="Gene expression", bty="n", main=main)
	abline(h=c(0,1), col="darkgray")
	polygon(x=c(12, 16, 16, 12), y=c(0,0,1,1), border=NA, col="lightgray")

	for (i in 1:nrow(W)) {
		lines(x=0:(ncol(mm$full)-1), mm$full[i,], col=col.genes[i], lwd=2)
		if (gene.names) 
			text(ncol(mm$full)-1.2, mm$full[i,ncol(mm$full)], rownames(W)[i], col=col.genes[i], pos=4)
	}
}

fit.plot <- function(thetas, s=rep(10, length(thetas)), phen=rep(NA, length(thetas)), xlab="Fitness", ylab="Gene expression", ...) {
	inv.gauss <- function(theta, s, reverse=FALSE, yrange = c(0,1)) {
	function(y) {
	ans <- if (reverse)
		theta + sqrt(-s*log(y))/s
	else
		theta - sqrt(-s*log(y))/s
	ans[ans < yrange[1]] <- NA
	ans[ans > yrange[2]] <- NA
	ans
	}}
	
	opar <- par(mar=c(par("mar")[1], 0, par("mar")[3], 2))
	
	plot(NULL, xlim=c(0,1.1), ylim=c(0,1), xlab=xlab, ylab=ylab, bty="n", xaxt="n", ...)
	axis(1, at=c(0,1/2,1), labels=c("0", "0.5", "1"))
	abline(h=c(0,1), col="darkgray")
	for (i in seq_along(thetas)) {
		if (is.na(thetas[i]) || is.na(s[i]) || s[i]==0)
			next

		curve(inv.gauss(thetas[i], s[i])(x), add=TRUE, col=col.genes[i], lwd=2)
		curve(inv.gauss(thetas[i], s[i], TRUE)(x), add=TRUE, col=col.genes[i], lwd=2)
		if (!is.na(phen[i]))
			lines(x=c(-0.3, 1), y=rep(phen[i],2), col=col.genes[i], lty=2, xpd=NA)
		points(rep(1, length(thetas)), thetas, pch=19, col="blue", cex=2)
		text(rep(1, length(thetas)), thetas, bquote(theta[.(names.genes[i])]), col="blue", pos=4, xpd=NA)
	}
	par(mar=opar)
}


thetas <- c(NA, 0.6, NA, NA)

W.before <- W.after <- rbind(
	c( 0, 0, 0, 0),
    c( 0.4, 0, -0.1, 0),
    c( 0, 0,  0, 0.5),
    c( 0, 0,  0, 0))
W.after[2,4] <- 0.6
rownames(W.before) <- rownames(W.after) <- colnames(W.before) <- colnames(W.after) <- names.genes
      

pdf("fig1.pdf", width=12, height=9)
layout(matrix(c(1:12), byrow=TRUE, ncol=4), widths=c(0.2,0.2,0.4,0.1))

par(cex=0.8, oma=c(0, 5, 0, 0))

W.plot(W.before)
mtext(side=2, text="Original\nNetwork", line=4, cex=2)
net.plot(W.before)
dyn.plot(W.before, env=env.default)
fit.plot(thetas, ylab="", yaxt="n", phen=model.M2(W.before, env=env.default)$mean)

W.plot(W.before, bold.caption=1)
mtext(side=2, text="Environment\nchange", line=4, cex=2, col="darkgreen")
net.plot(W.before, font.caption=c(2,1,1,1))
dyn.plot(W.before, env=env.alt)
fit.plot(thetas, ylab="", yaxt="n", phen=model.M2(W.before, env=env.alt)$mean)

W.plot(W.after, bold=which(W.before != W.after))
mtext(side=2, text="Mutant\nNetwork", line=4, cex=2, col="red")

net.plot(W.after, bold=which(W.before != W.after))
dyn.plot(W.after, env=env.default)
fit.plot(thetas, ylab="", yaxt="n", phen=model.M2(W.after, env=env.default)$mean)
dev.off()
