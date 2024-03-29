# Plotting helper function
# sourcing common-fig.R should be enough, other common functions are called from here.

source("./common-precalc.R")

source("../src/relvfit.R")
source("../src/reactnorm.R")
source("../src/connect.R")
source("../src/network.R")
source("../src/evolGW.R")
source("../src/variation.R")


# List of user-friendly functions
# (x: scalar, v: scalar or vector) 
# more options are available in most cases
#
### legname(v)        : 
#        transform simulation code into a user-friendly caption (e.g. "nosel" -> "Drift")
### generation.axis() :
#        add a properly formatted x axis for generations (don't forget to plot with xaxt="n")
### subpanel(x)       :
#        add a subpanel caption to the top-left corner
### avgcol(c1,c2)             : helper function, computes the average between two colors
### lighten.col(col,factor)   : helper function, lighten a color
### makeTransparent(col,alpha): helper function, makes a color transparent
#
# Data plots
#
### plot.N(x)         :
#        dynamics of census and effective population sizes for simulation x
#        used in fig7, figS4, figS5, figS9
### plot.fitness(x)   :
#        dynamics of average fitness for simulation x
#        used in fig2, figS4, figS5, figS9
### plot.var(v, what)  
### plot.var.pheno(v)       :
#        dynamics of expression variance for a set of simulations v
#        (in the figures, what is always "expression")
#        plot.var used in fig7, figS4, figS5, figS9
#        plot.var.pheno in fig4B
### plot.var.neutral(v):
#        dynamics of the molecular variance, focusing on quasi-neutral sites
#        used in fig4A, fig7, figS4, figS5, figS9
### plot.evol(v):
#        evolutionary change (evolution speed of the genotype)
#		 used in fig3A
### plot.norm(v)      :
#        dynamics of the absolute value of reaction norms for plastic genes for simulations v
#        used in fig3B; fig7, figS4, figS5, figS9
### plot.nconn(v)     :
#        dynamics of the number of connections for simulations v
#        used in figS7A
### plot.inout.change(x, x.ref):
#        plots the number of in vs. out connections before and after domestication.
#        if x.ref is provided, adds a reference simulation (e.g. no domestication)
#        used in figS7B
### plot.inout.gainloss(v, deltaG):
#        dynamics of the number of gained or lost connections every deltaG generations
#        used in fig6A, fig7, figS4, figS5, figS9
### plot.network.feature(v, what):
#        dynamics of some network properties (what) for a set of simulations v
#        available so far: what="clusters"
#        used in fig6B, figS4, figS5, figS9
### plot.numconn.groups(nc, nc.ref):
#        plots the number of connections between different groups of genes. Requires a complex 
#        object, nc, calculated by mean.numconn.groups. Can also accept correlation groups.
#        if a reference is provided, the plot features the delta with the ref network.
#        used in figS8
### plot.Gdiff(v, deltaG):
#        dynamics of the speed of evolution for the G matrix (difference in G matrix every deltaG generations) 
#        used in figS4, figS5, figS6A, figS9
### plot.Gcor(v, deltaG):
#        dynamics of the average genetic correlation
#        used in fig5A, fig7
### plot.Gmat(x, gen) :
### plot.Gmat.all(x, gen) :
### plot.Gmat.legend():
#        average genetic correlation matrix at generation gen
#        used in fig5B, figS6B
  

first.gen <- function(mysim) {
	all.gen <- as.numeric(names(Ndyn.all[[mysim]]))
	mxx <- max(all.gen)
	ff <- mxx - show.gen
	all.gen[which.min(abs(all.gen-ff))]
}

generation.axis <- function(show.bottleneck=FALSE, mysim="default", ...) {
	mxx <- max(as.numeric(names(Ndyn.all[[mysim]])))
	first.g <- first.gen(mysim)
	if (show.bottleneck) {
		bd <- unlist(bottleneck.detect(Ndyn.all[[mysim]]))
		toshow <- sort(c(first.g, mxx, bd))
	} else {
		toshow <- mxx - pretty(c(0, mxx - first.g), n=3)
		toshow <- toshow[toshow > first.g]
	}
	axis(1, at=toshow, labels = toshow - mxx, ...)
}

subpanel <- function(x, adj=0.025, col="black", line=-1, outer=FALSE) {
	title(adj=adj, main=x, cex.main=1.4, col.main=col, line=line, outer=outer)
}

avgcol <- function(c1, c2) {
	rgb(t(mapply(c1,c2,FUN=function(x,y) colorRamp(c(x,y))(0.5))), max=255)
}

lighten.col <- function(colors, factor=0.5) {
    sapply(colors, function(col) { cc <- col2rgb(col); rgb(t(cc + (255 - cc)*factor),max=255) })
}

makeTransparent<-function(someColor, alpha=100)
{
	# From https://stackoverflow.com/questions/8047668/transparent-equivalent-of-given-color
	# Author: Nick Sabbe
	# Licence : CC-attribution-SA from the conditions of the website
	newColor<-col2rgb(someColor)
	apply(newColor, 2, function(curcoldata){rgb(red=curcoldata[1], green=curcoldata[2], blue=curcoldata[3],alpha=alpha, maxColorValue=255)})
}

# Total population size and effective population size
plot.N <- function(mysim, ylab="Population size", xlab="Generations", ylim=c(0, max(Ndyn.all[[mysim]])), show.quantiles=FALSE, print.Ne=TRUE, ...) {
	vfit  <-  relvfit.mean.cache(outdir.all[[mysim]])
	q1vfit <- relvfit.quantile.cache(outdir.all[[mysim]], quant=quantiles[1])
	q2vfit <- relvfit.quantile.cache(outdir.all[[mysim]], quant=quantiles[2])
	gen <-  meansim.all[[mysim]][,"Gen"]
	
	plot(NULL, xlim=c(first.gen(mysim), max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...) 
	
	nn <- Ndyn.all[[mysim]][as.character(gen)]
	
	lines(	x=as.numeric(names(Ndyn.all[[mysim]])), 
			y=Ndyn.all[[mysim]])
	
	Ne <- nn/(1+4*vfit)
	lines(	x=gen, 
			y=Ne, 
			col="blue")
	if(show.quantiles) {
		polygon(c(gen, rev(gen)), c(nn, rev(nn))/(1+4*c(q1vfit, rev(q2vfit))), border=NA, col=makeTransparent("blue"))
	}
	if (print.Ne) {
		my.Ne <- Ne[which(names(Ne) == as.character(first.gen(mysim))):length(Ne)]
		Ne.ref <- min(my.Ne[1], my.Ne[length(my.Ne)])
		gen.bottle <- range(which(my.Ne < 0.90*Ne.ref))
		Ne.bottle <- 1/(mean(1/my.Ne[gen.bottle[1]:gen.bottle[2]]))
		text(x=mean(as.numeric(names(my.Ne[gen.bottle]))), y=min(nn), paste0("Ne =\n", round(Ne.bottle)), pos=if (diff(as.numeric(names(my.Ne[gen.bottle]))) < 1000) 4 else 3)
	}
}

plot.fitness <- function(mysims, show.quantiles=FALSE, ylab="Fitness", xlab="Generations", ylim=c(0,1), lty=NULL, ...) {
	gen <-  meansim.all[[mysims[1]]][,"Gen"]
	
	plot(NULL, xlab=xlab, ylab=ylab, xlim=c(first.gen(mysims[1]), max(gen)), ylim=ylim, ...)

	for (mysim in mysims) {
		mfit <- meansim.all[[mysim]][,"MFit"]
		mfit[mfit < ylim[1]] <- ylim[1] # Avoids problems on log scales
		lines(gen, mfit, lty=if(is.null(lty)) lty.sce[mysim] else lty)
		if (show.quantiles) {
			q1 <- quantile.sim.cache(outdir.all[[mysim]], quant=quantiles[1], mc.cores=mc.cores)[,"MFit"]
			q2 <- quantile.sim.cache(outdir.all[[mysim]], quant=quantiles[2], mc.cores=mc.cores)[,"MFit"]
			q1[q1 < ylim[1]] <- ylim[1]
			q2[q2 < ylim[1]] <- ylim[1]
			polygon(c(gen, rev(gen)), c(q1, rev(q2)), col=makeTransparent("black"), border=NA)
		}
	}
}

# Variances

plot.var <- function(mysims, what=c("molecular", "expression")[1], ylim=NULL, xlab="Generation", ylab=NULL, y.factor=1, show.quantiles=FALSE, ...) {
	if (what == "molecular") {
		var.data <- sapply(mysims, function(mysim) molec.variation(meansim.all[[mysim]])[,-1], USE.NAMES=TRUE, simplify=FALSE)
		if (is.null(ylab)) ylab <- "Molecular variance"
	} else if (what == "expression") {
		var.data <- sapply(mysims, function(mysim) pheno.variation(meansim.all[[mysim]])[,-1], USE.NAMES=TRUE, simplify=FALSE)
		if (is.null(ylab)) ylab <- "Expression variance"
	} else {
		stop("plot.var: what= ", what, " is not a valid option.")
	}

	gen <- as.numeric(meansim.all[[mysims[1]]][,"Gen"])
	
	if (is.null(ylim)) ylim <- c(0, y.factor*max(unlist(var.data)))

	plot(NULL, xlim=c(first.gen(mysims[1]), max(gen)), ylim=ylim, ylab=ylab, xlab=xlab, ...)
	
	for (mysim in mysims) {
		col <- if (mysim %in% names(col.sce)) col.sce[mysim] else "black"
		lty <- if (mysim %in% names(lty.sce)) lty.sce[mysim] else 1
		yy <- rowMeans(var.data[[mysim]]*y.factor)
		my.yy <- mov.avg(yy, gen, size=window.avg)
		lines(as.numeric(names(my.yy)), my.yy, lty=lty, col=col)
		
		if (show.quantiles) {
			yy.q1 <- apply(var.data[[mysim]]*y.factor, 1, quantile, prob=quantiles[1])
			yy.q2 <- apply(var.data[[mysim]]*y.factor, 1, quantile, prob=quantiles[2])
			my.yy.q1 <- mov.avg(yy.q1, gen, size=window.avg)
			my.yy.q2 <- mov.avg(yy.q2, gen, size=window.avg)
			xx <- as.numeric(names(my.yy))
			polygon(c(xx, rev(xx)), c(my.yy.q1, rev(my.yy.q2)), col=makeTransparent(col), border=NA)
		}
	}
}


plot.var.pheno <- function(mysims,  ylim=NULL, xlab="Generation", ylab="Expression variance", y.factor=1, show.quantiles=FALSE, ...) {
	for (mysim in mysims) {
		var.data <- pheno.variation.mean.cache(outdir.all[[mysim]], mc.cores=mc.cores)
		gen <- as.numeric(names(var.data))
		col <- if (mysim %in% names(col.sce)) col.sce[mysim] else "black"
		lty <- if (mysim %in% names(lty.sce)) lty.sce[mysim] else 1
		if (mysim == mysims[1]) { # not very clean
			if (is.null(ylim)) ylim <- c(0, y.factor*max(unlist(var.data)))
			plot(NULL, xlim=c(first.gen(mysim), max(gen)), ylim=ylim, ylab=ylab, xlab=xlab, ...)
		}
		lines(gen, y.factor*var.data, lty=lty, col=col)
		if (show.quantiles) {
			q1 <- pheno.variation.quantile.cache(outdir.all[[mysim]], quant=quantiles[1], mc.cores=mc.cores)
			q2 <- pheno.variation.quantile.cache(outdir.all[[mysim]], quant=quantiles[2], mc.cores=mc.cores)
			polygon(c(gen, rev(gen)), c(y.factor*q1, rev(y.factor*q2)), border=NA, col=makeTransparent(col))
		}
	}
}



plot.var.neutral <- function(mysims,  ylim=NULL, xlab="Generation", ylab="Molecular variance", expr.thresh=0.1, algorithm=c("lowexpr", "cleanW")[1], y.factor=1, show.quantiles=FALSE, ...) {
	
	for (mysim in mysims) {
		var.data <- molec.variation.neutral.all.mean.cache(outdir.all[[mysim]], expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores)
		gen <- as.numeric(names(var.data))
		col <- if (mysim %in% names(col.sce)) col.sce[mysim] else "black"
		lty <- if (mysim %in% names(lty.sce)) lty.sce[mysim] else 1
		
		if (mysim == mysims[1]) { # not very clean
			if (is.null(ylim)) ylim <- c(0, y.factor*max(unlist(var.data)))
			plot(NULL, xlim=c(first.gen(mysim), max(gen)), ylim=ylim, ylab=ylab, xlab=xlab, ...)
		}
		lines(gen, y.factor*var.data, lty=lty, col=col)
		if (show.quantiles) {
			q1 <- molec.variation.neutral.all.quantile.cache(outdir.all[[mysim]], quant=quantiles[1], expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores)
			q2 <- molec.variation.neutral.all.quantile.cache(outdir.all[[mysim]], quant=quantiles[2], expr.thresh=expr.thresh, algorithm=algorithm, mc.cores=mc.cores)
			polygon(c(gen, rev(gen)), c(y.factor*q1, rev(y.factor*q2)), border=NA, col=makeTransparent(col))
		}
	}
}



plot.evol <- function(mysims, ylim=NULL, xlab="Generation", ylab="Evolutionary change", ...) {
	
	gen <-  as.numeric(meansim.all[[mysims[1]]][,"Gen"]) # Just for the x scaling
	
	plot(NULL, xlim=c(first.gen(mysims[[1]]), max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	
	for (mysim in mysims) {
		ev.genes <- mean.Wdiff.dyn.cache(outdir.all[[mysim]], deltaG, mc.cores=mc.cores)[,-1]
		
		lines(as.numeric(rownames(ev.genes)), rowMeans(ev.genes), lty=lty.sce[mysim], col=col.sce[mysim])
	}
}


# Reaction norm, several simulations possible
plot.norm <- function(mysims, ylim=c(0, 1.2), xlab="Generation", ylab="|Reaction norm|", lty=NULL, show.quantiles=FALSE, ...) {
	gen <-  as.numeric(meansim.all[[mysims[1]]][,"Gen"])
	
	plot(NULL, xlim=c(first.gen(mysims[1]), max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...)
		
	for (mysim in mysims) {
		mean.norm <-  reaction.norm.mean.cache(outdir.all[[mysim]], FUN.to.apply=abs, sliding=TRUE, window.size=window.norm, mc.cores=mc.cores)

		yy.pp <- rowMeans(mean.norm[,selpattern.all[[mysim]]=="pp",drop=FALSE])
		yy.ps <- rowMeans(mean.norm[,selpattern.all[[mysim]]=="ps",drop=FALSE])
		yy.pn <- rowMeans(mean.norm[,selpattern.all[[mysim]]=="pn",drop=FALSE])
		
		xx <- as.numeric(rownames(mean.norm))
		
		lines(xx, yy.pp, col=col.sel["p"], lty=if(is.null(lty)) lty.sce[mysim] else lty)
		lines(xx, yy.ps, col=col.sel["s"], lty=if(is.null(lty)) lty.sce[mysim] else lty)
		lines(xx, yy.pn, col=col.sel["n"], lty=if(is.null(lty)) lty.sce[mysim] else lty)
		
		if (show.quantiles) {
				q1.norm <- reaction.norm.quantile.cache(outdir.all[[mysim]], quant=quantiles[1], FUN.to.apply=abs, sliding=TRUE, window.size=window.norm, mc.cores=mc.cores)
				q2.norm <- reaction.norm.quantile.cache(outdir.all[[mysim]], quant=quantiles[2], FUN.to.apply=abs, sliding=TRUE, window.size=window.norm, mc.cores=mc.cores)

				yy.pp <- c(rowMeans(q1.norm[,selpattern.all[[mysim]]=="pp",drop=FALSE]), rev(rowMeans(q2.norm[,selpattern.all[[mysim]]=="pp",drop=FALSE])))
				yy.ps <- c(rowMeans(q1.norm[,selpattern.all[[mysim]]=="ps",drop=FALSE]), rev(rowMeans(q2.norm[,selpattern.all[[mysim]]=="ps",drop=FALSE])))
				yy.pn <- c(rowMeans(q1.norm[,selpattern.all[[mysim]]=="pn",drop=FALSE]), rev(rowMeans(q2.norm[,selpattern.all[[mysim]]=="pn",drop=FALSE])))
				
				polygon(c(xx, rev(xx)), yy.pp, border=NA, col=makeTransparent(col.sel["p"]))
				polygon(c(xx, rev(xx)), yy.ps, border=NA, col=makeTransparent(col.sel["s"]))
				polygon(c(xx, rev(xx)), yy.pn, border=NA, col=makeTransparent(col.sel["n"]))
		}
	}
}

# Plots the dynamics of the number of connections. 
plot.nconn <- function(mysims, ylim=NULL, xlab="Generation", ylab="Number of connections", lty=NULL, col=NULL, show.quantiles=FALSE, ...) {
	nconn.all <- sapply(mysims, function(mysim) {
		mean.connect  <- number.connections.mean.cache(outdir.all[[mysim]],  mc.cores=mc.cores)
#~ 		mov.avg(mean.connect,  as.numeric(names(mean.connect)),  size=window.avg, min.gen=0)
	}, USE.NAMES=TRUE, simplify=FALSE)
	
	plot(NULL, 
		xlim=c(first.gen(mysims[1]), max(as.numeric(names(nconn.all[[1]])))), 
		ylim=if(is.null(ylim)) c(0, max(unlist(nconn.all))) else ylim, 
		xlab=xlab, ylab=ylab, ...)
	for (mysim in mysims) {
		mycol <- if(is.null(col)) col.sce[mysim] else col
		lines(as.numeric(names(nconn.all[[mysim]])),  nconn.all[[mysim]],  col=mycol, lty=if(is.null(lty)) lty.sce[mysim] else lty)
		if (show.quantiles) {
			q1 <- number.connections.quantile.cache(outdir.all[[mysim]], quant=quantiles[1], mc.cores=mc.cores)
			q2 <- number.connections.quantile.cache(outdir.all[[mysim]], quant=quantiles[2], mc.cores=mc.cores)
			xx <- as.numeric(names(q1))
			polygon(c(xx, rev(xx)), c(q1, rev(q2)), col=makeTransparent(mycol), border=NA)
		}
	}
}


plot.inout.change <- function(mysim, mysim.ref=NULL, regimes=c("s","p","n"), xlab="Connections to the gene", ylab="Connections from the gene", xlim=NULL, ylim=NULL, ...) {
	# assumes that the choice of mysim and mysim.ref make sense...
	mean.inout <- function(x)
		list(connect.in = rowMeans(sapply(x, "[[", "connect.in")), connect.out=rowMeans(sapply(x, "[[", "connect.out")))
	
	genselchange <- selectionchange.detect(meansim.all[[mysim]])
	
	before.dom      <- inout.gen.cache(outdir.all[[mysim]], gen=genselchange, mc.cores=mc.cores)
	mean.before.dom <- mean.inout (before.dom)
	after.dom       <- inout.gen.cache(outdir.all[[mysim]], gen=meansim.all[[mysim]][nrow(meansim.all[[mysim]]),"Gen"], mc.cores=mc.cores)
	mean.after.dom  <- mean.inout (after.dom)

	mean.after.ref <- NULL
	if (!is.null(mysim.ref)) {
		mean.after.ref <- mean.inout(inout.gen.cache(outdir.all[[mysim.ref]], gen=meansim.all[[mysim.ref]][nrow(meansim.all[[mysim.ref]]),"Gen"], mc.cores=mc.cores))
	}
	
	if (is.null(xlim)) 
		xlim <- c(0, max(c(mean.before.dom$connect.in, mean.after.dom$connect.in)))
	if (is.null(ylim))
		ylim <- c(0, max(c(mean.before.dom$connect.out, mean.after.dom$connect.out)))
	
	first.sel  <- substr(selpattern.all[[mysim]], 1, 1)
	second.sel <- substr(selpattern.all[[mysim]], 2, 2)
	both.sel   <- unique(selpattern.all[[mysim]])[-1] # removing the cc case
	
	plot(NULL, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, ...)
	
	points(
		x=sapply(regimes, function(nn) mean(mean.before.dom$connect.in[-1][first.sel==nn])),
		y=sapply(regimes, function(nn) mean(mean.before.dom$connect.out[-1][first.sel==nn])), 
		pch=19, col=col.sel[regimes], cex=2)
	x1 <- sapply(both.sel, function(nn) mean(mean.after.dom$connect.in[-1][selpattern.all[[mysim]]==nn]))
	y1 <- sapply(both.sel, function(nn) mean(mean.after.dom$connect.out[-1][selpattern.all[[mysim]]==nn]))
	points(x=x1, y=y1, pch=17, col=col.sel[substr(both.sel,2,2)], cex=2)
	arrows(
		x0=sapply(both.sel, function(nn) mean(mean.before.dom$connect.in[-1][substr(selpattern.all[[mysim]],1,1)==substr(nn,1,1)])),
		x1=x1,
		y0=sapply(both.sel, function(nn) mean(mean.before.dom$connect.out[-1][substr(selpattern.all[[mysim]],1,1)==substr(nn,1,1)])),
		y1=y1,
		col="darkgray", length=0.1)
	if(!is.null(mysim.ref))
		points(
			x=sapply(regimes, function(nn) mean(mean.after.ref$connect.in[-1][selpattern.all[[mysim.ref]]==nn])),
			y=sapply(regimes, function(nn) mean(mean.after.ref$connect.out[-1][selpattern.all[[mysim.ref]]==nn])), 
			pch=1, col=col.sel[regimes], cex=2)
}

plot.inout.gainloss <- function(mysims, deltaG=1, xlab="Generation", ylab="Nb connections", ylim=NULL, lty=NULL, show.quantiles=FALSE, ...) {
	
	iogl <- sapply(mysims, function(mysim) {
		fl <- delta.inout.mean.cache(outdir.all[[mysim]], deltaG, mc.cores=mc.cores)
	}, USE.NAMES=TRUE, simplify=FALSE)
	
	gen <-as.numeric(rownames(iogl[[1]]))
	if (is.null(ylim))
		ylim <- c(-1,1)*max(unlist(iogl))
			
	plot(NULL, xlim=c(first.gen(mysims[1]), max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	abline(h=0, col="darkgray", lty=3)
	
	for (mysim in mysims) {
		lines(gen, iogl[[mysim]][,"gain"], lty=if(is.null(lty)) lty.sce[mysim] else lty, col=col.gl["Gain"])
		lines(gen, -iogl[[mysim]][,"loss"], lty=if(is.null(lty)) lty.sce[mysim] else lty, col=col.gl["Loss"])
		if (show.quantiles) {
			q1 <- delta.inout.quantile.cache(outdir.all[[mysim]], quant=quantiles[1], deltaG=deltaG, mc.cores=mc.cores)
			q2 <- delta.inout.quantile.cache(outdir.all[[mysim]], quant=quantiles[2], deltaG=deltaG, mc.cores=mc.cores)
			xx <- as.numeric(rownames(q1))
			polygon(c(xx, rev(xx)), c(q1[,"gain"], rev(q2[,"gain"])), col=makeTransparent(col.gl["Gain"]), border=NA)
			polygon(c(xx, rev(xx)), c(-q1[,"loss"], rev(-q2[,"loss"])), col=makeTransparent(col.gl["Loss"]), border=NA)
		}
	}
}

plot.network.feature <- function(mysims, what=c("path"), directed=TRUE, deltaG=1, ylab=what, xlab="Generation", ylim=NULL, lty=NULL, col=NULL, show.quantiles=FALSE, ...) {
	callFUN <- if (what == "path") {
			'function(w) igraph::average.path.length(igraph::graph_from_adjacency_matrix(sign(cleanW(w))))'
		} else if (what == "diameter") {
			'function(w) igraph::diameter(igraph::graph_from_adjacency_matrix(sign(cleanW(w))))'
		} else if (what == "clusters") {
			'function(w) igraph::clusters(igraph::graph_from_adjacency_matrix(sign(cleanW(w))))$no'
		} else if (what == "centralization") {
			'function(w) igraph::centr_degree(igraph::graph_from_adjacency_matrix(sign(cleanW(w))))$centralization'
		} else if (what == "modularity") {
			'function(w) igraph::modularity(igraph::cluster_walktrap(igraph::graph_from_adjacency_matrix(sign(cleanW(w)))))'
		} else {
			stop()
		}

	netf.all <- sapply(mysims, function(mysim) {
		list(mean = WFUN.mean.cache(outdir.all[[mysim]], WFUN=callFUN, deltaG=deltaG, mc.cores=mc.cores, cache.subdir=paste0("Rcache-", what)),
					q1    = if (show.quantiles) WFUN.quantile.cache(outdir.all[[mysim]], quant=quantiles[1], WFUN=callFUN, deltaG=deltaG, mc.cores=mc.cores, cache.subdir=paste0("Rcache-", what)) else NA,
					q2    = if (show.quantiles) WFUN.quantile.cache(outdir.all[[mysim]], quant=quantiles[2], WFUN=callFUN, deltaG=deltaG, mc.cores=mc.cores, cache.subdir=paste0("Rcache-", what)) else NA)
	}, USE.NAMES=TRUE, simplify=FALSE)

	plot(NULL, 
		xlim=c(first.gen(mysims[1]), max(as.numeric(names(netf.all[[1]]$mean)))), 
		ylim=if(is.null(ylim)) c(0, max(unlist(netf.all), na.rm=TRUE)) else ylim, 
		xlab=xlab, ylab=ylab, ...)
	
	for (mysim in mysims) {
		xx <- as.numeric(names(netf.all[[mysim]]$mean))
		mycol <- if(is.null(col)) col.sce[mysim] else col
		lines(xx, netf.all[[mysim]]$mean,  col=mycol, lty=if(is.null(lty)) lty.sce[mysim] else lty)
		if (show.quantiles) {
			polygon(c(xx, rev(xx)), c(netf.all[[mysim]]$q1, rev(netf.all[[mysim]]$q2)), col=makeTransparent(mycol), border=NA)
		}
	}
}

plot.numconn.groups <- function(numconn, group.names=colnames(numconn$plus), 
								numconn.ref=NULL, directed=TRUE, ann.text=TRUE, col.scale.plus=NULL, col.scale.minus=NULL, 
								lwd.arr=2, pos.shift.plus=0.80, pos.shift.minus=0.60, thresh=0.05, ...) {
													
	circ.arc <- function(theta1=0, theta2=2*pi, n=100) { tt <- seq(theta1, theta2, length.out=n); cbind(cos(tt), sin(tt)) }
	posit.angle <- function(angle) { angle <- angle %% (2*pi); if (angle > pi/2 && angle <= 3*pi/2) angle <- angle + pi; angle %% (2*pi)}
	if (is.null(col.scale.plus))
		col.scale.plus <- colorRampPalette(c("white","black"))(100)
	if (is.null(col.scale.minus))
		col.scale.minus <- colorRampPalette(c("white","red"))(100)
	
	delta.angle <- 0.25 # angle between two arrows
	arr.dist  <- 0.15 # distance between the group name and the arrows
	self.angle <- 1.4*pi
	ann.text.options <- list(
		pos.shift.plus=pos.shift.plus, 
		pos.shift.minus=pos.shift.minus, 
		text.cex=0.7, 
		col.plus=rev(col.scale.plus)[1], 
		col.minus=rev(col.scale.minus)[1], 
		thresh=thresh, 
		digits=2)
	
	par(mar=c(0.1,0.1,2,0.1))
	plot(NULL, xlim=c(-1.2,1.2), ylim=c(-1.2,1.2), axes=FALSE, ann=FALSE, asp=1, ...)
	lg <- length(group.names)
	xy.groups <- cbind(cos(2*pi/lg*(0:(lg-1))), sin(2*pi/lg*(0:(lg-1))))
	
	if (is.null(names(group.names))) names(group.names) <- group.names
	# plots the groups
	text(xy.groups[,1], xy.groups[,2], parse(text=names(group.names)))
	
	numconn$plus[!is.finite(numconn$plus)] <- 0
	numconn$minus[!is.finite(numconn$minus)] <- 0
	
	for (i in 1:lg) {
		for (j in 1:lg) {
			if (i != j) {
				# angle between i and j
				alpha <- atan((xy.groups[j,2]-xy.groups[i,2])/(xy.groups[j,1]-xy.groups[i,1]))
				if (xy.groups[i,1] > xy.groups[j,1]) alpha <- alpha - pi
				
				x.i.plus  <- xy.groups[i,1]+arr.dist*cos(alpha-delta.angle/2)
				y.i.plus  <- xy.groups[i,2]+arr.dist*sin(alpha-delta.angle/2)
				x.i.minus <- xy.groups[i,1]+arr.dist*cos(alpha-3*delta.angle/2)
				y.i.minus <- xy.groups[i,2]+arr.dist*sin(alpha-3*delta.angle/2)
				
				x.j.plus  <- xy.groups[j,1]+arr.dist*cos(pi+alpha+delta.angle/2)
				y.j.plus  <- xy.groups[j,2]+arr.dist*sin(pi+alpha+delta.angle/2)
				x.j.minus <- xy.groups[j,1]+arr.dist*cos(pi+alpha+3*delta.angle/2)
				y.j.minus <- xy.groups[j,2]+arr.dist*sin(pi+alpha+3*delta.angle/2)
				
				col.plus  <- col.scale.plus [round(numconn$plus[j,i]* length(col.scale.plus))]
				col.minus <- col.scale.minus[round(numconn$minus[j,i]*length(col.scale.minus))]
				
				if (directed) {
					arrows(x0=x.i.plus,  x1=x.j.plus,  y0=y.i.plus,  y1=y.j.plus,  length=0.1,  col=col.plus,  lwd=lwd.arr)
					arrows(x0=x.i.minus, x1=x.j.minus, y0=y.i.minus, y1=y.j.minus, length=0.05, col=col.minus, lwd=lwd.arr, angle=90)
				} else { # The non-directed case is a bit of a hack, reusing existing variables
					if (i > j) { # otherwise no need to draw anything
						arrows(x0=x.i.plus,  x1=x.j.plus,  y0=y.i.plus,  y1=y.j.plus, code=3, length=0,  col=col.plus,  lwd=lwd.arr)
					} else { # j < i necessarily, as the case i==j is treated elsewhere
						arrows(x0=x.i.plus,  x1=x.j.plus,  y0=y.i.plus,  y1=y.j.plus, code=3, length=0,  col=col.minus,  lwd=lwd.arr)
					}
				}
				
				if (ann.text) {
#~ 					if (numconn$plus[j,i] > ann.text.options$thresh) {
						txt.num <- if (!directed && j > i) {
							if (is.null(numconn.ref))
								numconn$minus[i,j]
							else
								numconn$minus[i,j] - numconn.ref$minus[i,j]
						} else {
							if (is.null(numconn.ref))
								numconn$plus[j,i]
							else
								numconn$plus[j,i] - numconn.ref$plus[j,i]
						}
						txt.str <- sprintf(paste0("%", if (is.null(numconn.ref)) "" else "+", "1.", ann.text.options$digits, "f"), txt.num)
						if (abs(txt.num) > ann.text.options$thresh)
						text(x=(1-ann.text.options$pos.shift.plus)*x.i.plus+ann.text.options$pos.shift.plus*x.j.plus, 
							y=(1-ann.text.options$pos.shift.plus)*y.i.plus+ann.text.options$pos.shift.plus*y.j.plus, 
							txt.str, 
							cex=ann.text.options$text.cex, 
							col=if(!directed && j > i) ann.text.options$col.minus else ann.text.options$col.plus,
							srt=180*(posit.angle(alpha)/pi))
					}
#~ 					if (numconn$minus[j,i] > ann.text.options$thresh && directed) {
					if (directed) {
						txt.num <- if (is.null(numconn.ref)) numconn$minus[j,i]	else numconn$minus[j,i] - numconn.ref$minus[j,i]
						txt.str <- sprintf(paste0("%", if (is.null(numconn.ref)) "" else "+", "1.", ann.text.options$digits, "f"), txt.num)
						if (abs(txt.num) > ann.text.options$thresh)
						text(x=(1-ann.text.options$pos.shift.minus)*x.i.minus+ann.text.options$pos.shift.minus*x.j.minus, 
							y=(1-ann.text.options$pos.shift.minus)*y.i.minus+ann.text.options$pos.shift.minus*y.j.minus, 
							txt.str, 
							cex=ann.text.options$text.cex, 
							col=ann.text.options$col.minus,
							srt=180*(posit.angle(alpha)/pi))
					}
#~ 				}
			} else { # i == j
				# angle of i around the circle
				alpha <- (i-1)*2*pi/lg
				# the regulation circle is opposite to the center
				cc <- circ.arc(alpha-self.angle/2, alpha+self.angle/2)
				cc.plus <- t(t(cc)*arr.dist+(1+0.5*arr.dist)*xy.groups[i,])
				cc.minus <- t(t(cc)*(1-delta.angle)*arr.dist+(1+0.5*arr.dist)*xy.groups[i,])
				col.plus <- col.scale.plus[round(numconn$plus[i,i]*length(col.scale.plus))]
				col.minus <- col.scale.minus[round(numconn$minus[i,i]*length(col.scale.minus))]
				lines(cc.plus[1:(nrow(cc.plus)-1), 1], cc.plus[1:(nrow(cc.plus)-1),2], col=col.plus, lty=1, lwd=lwd.arr)
				if (directed)
					arrows(x0=cc.plus[nrow(cc.plus)-1,1], x1=cc.plus[nrow(cc.plus),1], y0=cc.plus[nrow(cc.plus)-1,2], y1=cc.plus[nrow(cc.plus),2], col=col.plus, length=0.1, lwd=lwd.arr)
				lines(cc.minus[2:nrow(cc.minus), 1], cc.minus[2:nrow(cc.minus),2], col=col.minus, lty=1, lwd=lwd.arr)
				if (directed)
					arrows(x0=cc.minus[2,1], x1=cc.minus[1,1], y0=cc.minus[2,2], y1=cc.minus[1,2], col=col.minus, length=0.05, angle=90, lwd=lwd.arr)
				
				if (ann.text) {
#~ 					if (numconn$plus[i,i] > ann.text.options$thresh) {
						txt.num <- if (is.null(numconn.ref)) numconn$plus[i,i]	else numconn$plus[i,i] - numconn.ref$plus[i,i]
						txt.str <- sprintf(paste0("%", if (is.null(numconn.ref)) "" else "+", "1.", ann.text.options$digits, "f"), txt.num)
						if (is.finite(txt.num) && abs(txt.num) > ann.text.options$thresh)
						text(x=cc.plus[round(ann.text.options$pos.shift.plus*nrow(cc.plus)),1], 
							y=cc.plus[round(ann.text.options$pos.shift.plus*nrow(cc.plus)),2], 
							txt.str, 
							cex=ann.text.options$text.cex, 
							col=ann.text.options$col.plus)
#~ 					}
#~ 					if (numconn$minus[i,i] > ann.text.options$thresh) {
						txt.num <- if (is.null(numconn.ref)) numconn$minus[i,i]	else numconn$minus[i,i] - numconn.ref$minus[i,i]
						txt.str <- sprintf(paste0("%", if (is.null(numconn.ref)) "" else "+", "1.", ann.text.options$digits, "f"), txt.num)						
						if (is.finite(txt.num) && abs(txt.num) > ann.text.options$thresh)
						text(x=cc.minus[round(ann.text.options$pos.shift.minus*nrow(cc.minus)),1], 
							y=cc.minus[round(ann.text.options$pos.shift.minus*nrow(cc.minus)),2], 
							txt.str, 
							cex=ann.text.options$text.cex, 
							col=ann.text.options$col.minus)
#~ 					}
				}
			}
		}
	}
}


plot.Gdiff <- function(mysims, show.quantiles=FALSE, deltaG=1, xlab="Generation", ylab="Change in G matrix", ylim=c(0,1), col=NULL, lty=NULL, ...) {
	gd <- sapply(mysims, function(mysim) {
		list(
			mean = Gdiff.mean.cache(outdir.all[[mysim]], deltaG, mc.cores=mc.cores),
			q1   = if (show.quantiles) Gdiff.quantile.cache(outdir.all[[mysim]], quant=quantiles[1], deltaG, mc.cores=mc.cores) else NA,
			q2   = if (show.quantiles) Gdiff.quantile.cache(outdir.all[[mysim]], quant=quantiles[2], deltaG, mc.cores=mc.cores) else NA)
	}, USE.NAMES=TRUE, simplify=FALSE)
	
	gen <-as.numeric(names(gd[[1]]$mean))

	plot(NULL, xlim=c(first.gen(mysims[1]), max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	
	for (mysim in mysims) {
		mycol <- if(is.null(col)) col.sce[mysim] else col
		lines(gen, gd[[mysim]]$mean, lty=if(is.null(lty)) lty.sce[mysim] else lty, col=mycol)
		if (show.quantiles) {
			polygon(c(gen, rev(gen)), c(gd[[mysim]]$q1, rev(gd[[mysim]]$q2)), col=makeTransparent(mycol),  border=NA)
		}
	}
}

plot.Gcor <- function(mysims, xlab="Generations", ylab="Genetic correlations", ylim=NULL, show.quantiles=FALSE, ...) {
	
	gen <- as.numeric(names(Gcor.mean.cache(outdir.all[[mysims[1]]], mc.cores=mc.cores)))
	
	plot(NULL, xlim=c(first.gen(mysims[1]), max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	
	for (mysim in mysims) {
		gc <- Gcor.mean.cache(outdir.all[[mysim]], mc.cores=mc.cores)
		col <- if (mysim %in% names(col.sce)) col.sce[mysim] else "black"
		lty <- if (mysim %in% names(lty.sce)) lty.sce[mysim] else 1
		lines(gen, gc, lty=lty, col=col)
		if (show.quantiles) {
			q1 <- Gcor.quantile.cache(outdir.all[[mysim]], quant=quantiles[1], mc.cores=mc.cores)
			q2 <- Gcor.quantile.cache(outdir.all[[mysim]], quant=quantiles[2], mc.cores=mc.cores)
			
			polygon(c(gen, rev(gen)), c(q1, rev(q2)), border=NA, col=makeTransparent(col))
		}
	}
}

plot.Gmat <- function(mysim, gen, absolute=TRUE, cols=NULL, ...) {
	.transfG <- function(G)  if (absolute) abs(cov2cor(G)) else cov2cor(G)
	
	stopifnot(length(gen) <= 2) # if one gen: displays a cor matrix, if 2 gens: below and above diag. 
	
	if(is.null(cols))
		cols <- colorRampPalette(col.cor[(if(absolute) 2 else 1):length(col.cor)])(1001)

	st <- if (absolute) 0 else -1
	
	Glist1 <- Ggen.files.cache(outdir.all[[mysim]], gen[1], mc.cores=mc.cores)
	Glist1 <- lapply(Glist1, .transfG)
	
	if (length(gen) == 2) {
		Glist2 <- Ggen.files.cache(outdir.all[[mysim]], gen[2], mc.cores=mc.cores)
		Glist2 <- lapply(Glist2, .transfG)
	}
	
	Gmean <- rowMeans(do.call(abind, c(Glist1, list(along=3))), dims=2)
	if (length(gen) == 2) {
		Gmean2 <- rowMeans(do.call(abind, c(Glist2, list(along=3))), dims=2)
		Gmean[upper.tri(Gmean)] <- Gmean2[upper.tri(Gmean2)]
	}
	
	ng <- ncol(Gmean)-1
	
	image(x=1:ng, y=1:ng, t(Gmean[(ng+1):2,2:(ng+1)]), axes=FALSE, xlab="", ylab="", zlim=c(st, 1), col=cols, ...)
}

plot.Gmat.legend <- function(absolute=TRUE, ylab="Genetic correlation", cols=NULL) {
	if(is.null(cols))
		cols <- colorRampPalette(col.cor[(if(absolute) 2 else 1):length(col.cor)])(1001)
	
	st <- if (absolute) 0 else -1
	yl <- paste0(if (absolute) "|" else "", ylab, if (absolute) "|" else "")
	image(y=seq(st, 1, length.out=length(cols)), z=t(as.matrix(seq(st, 1, length.out=length(cols)-1))), zlim=c(st,1), col=cols, xaxt="n", xlab="", ylab=yl)
}

plot.Gmat.all <- function(mysim, gens, absolute=TRUE, cols, sel1, sel2=sel1, pch=c(same=1, diff=2), ylim=NULL, xlab="", ylab=NA, ...) {
	.transfG <- function(G)  if (absolute) abs(cov2cor(G)) else cov2cor(G)
		
	stopifnot(length(gens) == 2)

	st <- if (absolute) 0 else -1
	if (is.na(ylab)) ylab <- paste0(if(absolute)  "|", "Genetic correlation", if(absolute) "|")
		
	Glist1 <- Ggen.files.cache(outdir.all[[mysim]], gens[1], mc.cores=mc.cores)
	Glist1 <- lapply(Glist1, .transfG)
	Gmean1 <- rowMeans(do.call(abind, c(Glist1, list(along=3))), dims=2)[-1,-1]
	
	Glist2 <- Ggen.files.cache(outdir.all[[mysim]], gens[2], mc.cores=mc.cores)
	Glist2 <- lapply(Glist2, .transfG)
	Gmean2 <- rowMeans(do.call(abind, c(Glist2, list(along=3))), dims=2)[-1,-1]

	ng <- ncol(Gmean1)

	plot(NULL, xlim=c(-0.2,1.2), ylim=if(is.null(ylim)) c(st, 1) else ylim, xaxt="n", xlab=xlab, ylab=ylab, ...)
	
	for (i1 in 1:(ng-1))
		for (i2 in (i1+1):ng) {
			avgcol1 <- avgcol(cols[sel1[i1]], cols[sel1[i2]])
			avgcol2 <- avgcol(cols[sel2[i1]], cols[sel2[i2]])
			dx <- rnorm(1,0,sd=0.02)
			lines(x=dx+c(0,1), y=c(Gmean1[i1,i2], Gmean2[i1,i2]), lty=1, col=avgcol(avgcol1, avgcol2))
			points(x=dx+0, y=Gmean1[i1,i2], col=avgcol1, pch=if(sel1[i1]==sel1[i2]) pch[1] else pch[2])
			points(x=dx+1, y=Gmean2[i1,i2], col=avgcol2, pch=if(sel2[i1]==sel2[i2]) pch[1] else pch[2])	
		}
	axis(1, at=c(0,1), c("Before domestication", "Present"))
}
