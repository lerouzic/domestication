# Plotting helper function
# sourcing common-fig.R should be enough, other common functions are called from here.

source("./common-precalc.R")

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
### plot.N(x)         :
#        dynamics of census and effective population sizes for simulation x
### plot.fitness(x)   :
#        dynamics of average fitness for simulation x
### plot.var(v, what)       :
#        dynamics of variance for a set of simulations v
#        what can be "molecular" or "expression"
### plot.var.gene(x, what)  :
#        dynamics of variance for simulation x, one color per gene group
#        what can be "molecular" or "expression"
### plot.norm(v)      :
#        dynamics of the absolute value of reaction norms for plastic genes for simulations v
### plot.nconn(v)     :
#        dynamics of the number of connections for simulations v
### plot.inout.change(x, x.ref):
#        plots the number of in vs. out connections before and after domestication.
#        if x.ref is provided, adds a reference simulation (e.g. no domestication)
### plot.inout.gainloss(v, deltaG):
#        dynamics of the number of gained or lost connections every deltaG generations
### plot.network.feature(v, what):
#        dynamics of some network properties (what) for a set of simulations v
#        available so far: what="nbconn" | "modularity"
### plot.numconn.groups(nc, nc.ref):
#        plots the number of connections between different groups of genes. Requires a complex 
#        object, nc, calculated by mean.numconn.groups. Can also accept correlation groups.
#        if a reference is provided, the plot features the delta with the ref network.
### plot.Gdiff(v, deltaG):
#        dynamics of the speed of evolution for the G matrix (difference in G matrix every deltaG generations) 
### plot.GPC(v)       :
#        dynamics of the part of the genetic variance explained by the first PC
### plot.Gmat(x, gen) :
#        average genetic correlation matrix at generation gen
### plot.Gmat.legend():
#        color scale corresponding to plot.Gmat
#      


# Consistent caption text
legname <- function(nn) {
	legn <- c(
		default ="Default",
		nobot   ="No bottleneck",
		noselc  ="Constant selection",
		nosel   ="Drift",
		nomut   ="No new mutations",
		smallsel="Less selected genes",
		strongsel="Strong selection",
		strongbot="Strong bottleneck",
		largenet="Large network")
		
	ifelse(nn %in% names(legn), legn[nn], nn)
}

generation.axis <- function(show.bottleneck=FALSE, ...) {
	mxx <- max(as.numeric(names(Ndyn.all[["default"]])))
	if (show.bottleneck) {
		bd <- unlist(bottleneck.detect(Ndyn.all[["default"]]))
		toshow <- sort(c(first.gen, mxx, bd))
	} else {
		toshow <- mxx - pretty(c(0, mxx - first.gen), n=3)
		toshow <- toshow[toshow > first.gen]
	}
	axis(1, at=toshow, labels = toshow - mxx, ...)
}

subpanel <- function(x) {
	title(outer=FALSE,adj=0.025 ,main=x,cex.main=1.4,col="black",line=-1)
}


# Total population size and effective population size
plot.N <- function(mysim, ylab="Population size", xlab="Generations", ylim=c(0, max(Ndyn.all[[mysim]])), ...) {
	mfit <- meansim.all[[mysim]][,"MFit"]
	vfit <- meansim.all[[mysim]][,"VFit"]
	gen <-  meansim.all[[mysim]][,"Gen"]
	
	plot(NULL, xlim=c(first.gen, max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...) 
	
	lines(	x=as.numeric(names(Ndyn.all[[mysim]])), 
			y=Ndyn.all[[mysim]])
			
	lines(	x=gen, 
			y=Ndyn.all[[mysim]][as.character(gen)]/(1+4*vfit/mfit/mfit), 
			col="blue")
}

plot.fitness <- function(mysim, ylab="Fitness", xlab="Generations", ylim=c(0,1), ...) {
	mfit <- meansim.all[[mysim]][,"MFit"]
	gen <-  meansim.all[[mysim]][,"Gen"]
		
	plot(NULL, xlab=xlab, ylab=ylab, xlim=c(first.gen, max(gen)), ylim=ylim, ...)
	lines(gen, mfit)
}

# Variances

plot.var <- function(mysims, what=c("molecular", "expression")[1], ylim=NULL, xlab="Generation", ylab=NULL, y.factor=1, ...) {
	if (what == "molecular") {
		var.data <- sapply(mysims, function(mysim) molec.variation(meansim.all[[mysim]])[,-1], USE.NAMES=TRUE, simplify=FALSE)
		if (is.null(ylab)) ylab <- "Molecular variance"
	} else if (what == "expression") {
		var.data <- sapply(mysims, function(mysim) pheno.variation(meansim.all[[mysim]])[,-1], USE.NAMES=TRUE, simplify=FALSE)
		if (is.null(ylab)) ylab <- "Expression variance"
	} else {
		stop("plot.var: what= ", what, " is not a valid option.")
	}

	gen <- as.numeric(meansim.all[[mysims[[1]]]][,"Gen"])
	
	if (is.null(ylim)) ylim <- c(0, y.factor*max(unlist(var.data)))

	plot(NULL, xlim=c(first.gen, max(gen)), ylim=ylim, ylab=ylab, xlab=xlab, ...)
	
	for (mysim in mysims) {
		yy <- rowMeans(var.data[[mysim]]*y.factor)
		my.yy <- mov.avg(yy, gen, size=window.avg)
		lines(as.numeric(names(my.yy)), my.yy, lty=lty.sce[mysim], col=col.sce[mysim])
	}
}

plot.var.gene <- function(mysim, what=c("molecular", "expression")[1], ylim=NULL, xlab="Generation", ylab=NULL, y.factor=1, ...) {
	if (what == "molecular") {
		var.data <- molec.variation(meansim.all[[mysim]])[,-1]
		if (is.null(ylab)) ylab <- "Molecular variance"
	} else if (what == "expression") {
		var.data <- pheno.variation(meansim.all[[mysim]])[,-1]
		if (is.null(ylab)) ylab <- "Expression variance"
	} else {
		stop("plot.var: what= ", what, " is not a valid option.")
	}
	
	gen <- as.numeric(meansim.all[[mysim]][,"Gen"])

	sel.change.gen <- try(selectionchange.detect(meansim.all[[mysim]]), silent=TRUE)
	if (class(sel.change.gen) == "try-error") sel.change.gen <- max(gen)
	sel.before <- substr(selpattern.all[[mysim]], 1, 1)
	sel.before[sel.before == "c"] <- "s" # No need to distinguish constant and stable? 

	if (is.null(ylim)) ylim <- c(0, y.factor*max(var.data))

	plot(NULL, xlim=c(first.gen, max(gen)), ylim=ylim, ylab=ylab, xlab=xlab, ...)
	
	for (cc in unique(sel.before)) {
		yy <- (rowMeans(var.data[,sel.before==cc])*y.factor)[gen <= sel.change.gen]
		my.yy <- mov.avg(yy, gen[gen <= sel.change.gen], size=window.avg, min.gen=0)
		lines(as.numeric(names(my.yy)), my.yy, lty=1, col=col.sel[cc])
	}
	if (sel.change.gen < max(gen))
		for (cc in unique(selpattern.all[[mysim]])) {
			yy <- (rowMeans(var.data[,selpattern.all[[mysim]]==cc,drop=FALSE])*y.factor)[gen > sel.change.gen]
			my.yy <- mov.avg(yy, gen[gen > sel.change.gen], size=window.avg, min.gen=0)
			lines(as.numeric(names(my.yy)), my.yy, lty=lty.sel[substr(cc,1,1)], col=col.sel[substr(cc,2,2)])
		}
}


# Reaction norm, several simulations possible
plot.norm <- function(mysims, ylim=c(0, 1.2), xlab="Generation", ylab="|Reaction norm|", ...) {
	gen <-  as.numeric(meansim.all[[mysims[1]]][,"Gen"])
	
	plot(NULL, xlim=c(first.gen, max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...)
		
	for (mysim in mysims) {
		mean.norm <-  mean.norm.cache(outdir.all[[mysim]], FUN.to.apply=abs, sliding=TRUE, window.size=window.norm, mc.cores=mc.cores)

		yy.pp <- rowMeans(mean.norm[,selpattern.all[[mysim]]=="pp"])
		yy.ps <- rowMeans(mean.norm[,selpattern.all[[mysim]]=="ps"])
		yy.pn <- rowMeans(mean.norm[,selpattern.all[[mysim]]=="pn"])
		
		lines(as.numeric(rownames(mean.norm)), yy.pp, col=col.sel["p"], lty=lty.sce[mysim])
		lines(as.numeric(rownames(mean.norm)), yy.ps, col=col.sel["s"], lty=lty.sce[mysim])
		lines(as.numeric(rownames(mean.norm)), yy.pn, col=col.sel["n"], lty=lty.sce[mysim])
	}
}

# Plots the dynamics of the number of connections. 
plot.nconn <- function(mysims, ylim=NULL, xlab="Generation", ylab="Nb connections", ...) {
	nconn.all <- sapply(mysims, function(mysim) {
		mean.connect  <- mean.connect.cache(outdir.all[[mysim]],  mc.cores=mc.cores)
		mov.avg(mean.connect,  as.numeric(names(mean.connect)),  size=window.avg, min.gen=0)
	}, USE.NAMES=TRUE, simplify=FALSE)
	
	plot(NULL, 
		xlim=c(first.gen, max(as.numeric(names(nconn.all[[1]])))), 
		ylim=c(0, max(unlist(nconn.all))), 
		xlab=xlab, ylab=ylab, ...)
	for (mysim in mysims) {
		lines(as.numeric(names(nconn.all[[mysim]])),  nconn.all[[mysim]],  col=col.sce[mysim],  lty=lty.sce[mysim])
	}
}


plot.inout.change <- function(mysim, mysim.ref=NULL, regimes=c("s","p","n"), xlab="Connections to the gene", ylab="Connections from the gene", xlim=NULL, ylim=NULL, ...) {
	# assumes that the choice of mysim and mysim.ref make sense...
	mean.inout <- function(x)
		list(connect.in = rowMeans(sapply(x, "[[", "connect.in")), connect.out=rowMeans(sapply(x, "[[", "connect.out")))
	
	out.files.mysim <- list.files(pattern="out.*", path=list.dirs(outdir.all[[mysim]], full.names=TRUE, recursive=FALSE), full.names=TRUE)
	genselchange <- selectionchange.detect(meansim.all[[mysim]])
	
	before.dom      <- inout.gen.cache(out.files.mysim, gen=genselchange, epsilon=connect.threshold, env=connect.env, mc.cores=mc.cores)
	mean.before.dom <- mean.inout (before.dom)
	after.dom       <- inout.gen.cache(out.files.mysim, gen=meansim.all[[mysim]][nrow(meansim.all[[mysim]]),"Gen"], epsilon=connect.threshold, env=connect.env, mc.cores=mc.cores)
	mean.after.dom  <- mean.inout (after.dom)

	mean.after.ref <- NULL
	if (!is.null(mysim.ref)) {
		out.files.ref <- list.files(pattern="out.*", path=list.dirs(outdir.all[[mysim.ref]], full.names=TRUE, recursive=FALSE), full.names=TRUE)
		mean.after.ref <- mean.inout(inout.gen.cache(out.files.ref, gen=meansim.all[[mysim.ref]][nrow(meansim.all[[mysim.ref]]),"Gen"]))
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
	points(
		x=sapply(regimes, function(nn) mean(mean.after.ref$connect.in[-1][selpattern.all[[mysim.ref]]==nn])),
		y=sapply(regimes, function(nn) mean(mean.after.ref$connect.out[-1][selpattern.all[[mysim.ref]]==nn])), 
		pch=1, col=col.sel[regimes], cex=2)
}

plot.inout.gainloss <- function(mysims, deltaG=NA, xlab="Generation", ylab="Nb connections", ylim=NULL, ...) {
	
	iogl <- sapply(mysims, function(mysim) {
		out.files  <- list.files(pattern="out.*", path=list.dirs(outdir.all[[mysim]],  full.names=TRUE, recursive=FALSE), full.names=TRUE)
		fl <- mean.delta.inout.dyn.cache(out.files, deltaG)
	}, USE.NAMES=TRUE, simplify=FALSE)
	
	gen <-as.numeric(rownames(iogl[[1]]))
	if (is.null(ylim))
		ylim <- c(0, max(unlist(iogl)))
			
	plot(NULL, xlim=c(first.gen, max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	
	for (mysim in mysims) {
		lines(gen, iogl[[mysim]][,"gain"], lty=lty.sce[mysim], col=col.gl["gain"])
		lines(gen, iogl[[mysim]][,"loss"], lty=lty.sce[mysim], col=col.gl["loss"])		
	}	
}

plot.network.feature <- function(mysims, what=c("nbconn", "modularity")[1], algos=names(col.algo), ylab=NULL, xlab="Generation", ylim=NULL, ...) {
	library(igraph)
	comm.files <- function(out.files) {
	 mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		cc <- communities.dyn.cache(tt, epsilon=connect.threshold, env=env, directed=directed, mc.cores=max(1, round(mc.cores/4)))
		return(cc)
		}, mc.cores=min(4, mc.cores))
	}
	
	modn <- sapply(mysims, function(mysim) {
		out.files <- list.files(pattern="out.*", path=list.dirs(outdir.all[[mysim]], full.names=TRUE, recursive=FALSE), full.names=TRUE)
		comm      <- comm.files(out.files)
		sapply(algos, function(algo) 
			colMeans(do.call(rbind, mclapply(comm, function(cc) sapply(cc, function(ccc) {
				if (what=="nbconn")
					length(igraph::communities(ccc[[algo]]))
				else if (what=="modularity") 
					igraph::modularity(ccc[[algo]])
				else
					NA
			}), mc.cores=mc.cores))), simplify=FALSE)
		}, USE.NAMES=TRUE, simplify=FALSE)
		
	gen <- meansim.all[[mysims[1]]][,"Gen"]
	if (is.null(ylim)) 
		ylim <- c(0, max(unlist(modn)))
	if (is.null(ylab)) {
		if (what=="nbconn")     ylab <- "Number connections"
		if (what=="modularity") ylab <- "Modularity"
	}
	plot(NULL, xlim=c(first.gen, max(gen)), ylim=ylim, xlab="Generation", ylab=ylab, xaxt="n")
	for (cc in algos) {
		for (mysim in mysims) {
			my.modn  <- mov.avg(modn[[mysim]][[cc]],  as.numeric(names((modn[[mysim]][[cc]]))),  size=window.avg)
			lines(as.numeric(names(my.modn)),  my.modn,  col=col.algo[cc], lty=lty.sce[mysim])
		}
	}	
}

plot.numconn.groups <- function(numconn, group.names=colnames(numconn$plus), 
								numconn.ref=NULL, directed=TRUE, ann.text=TRUE, col.scale.plus=NULL, col.scale.minus=NULL, 
								lwd.arr=2, pos.shift.plus=0.80, pos.shift.minus=0.60, ...) {
													
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
		thresh=0.05, 
		digits=2)
	
	par(mar=c(0.1,0.1,4,0.1))
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
					if (numconn$plus[j,i] > ann.text.options$thresh) {
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
						text(x=(1-ann.text.options$pos.shift.plus)*x.i.plus+ann.text.options$pos.shift.plus*x.j.plus, 
							y=(1-ann.text.options$pos.shift.plus)*y.i.plus+ann.text.options$pos.shift.plus*y.j.plus, 
							txt.str, 
							cex=ann.text.options$text.cex, 
							col=if(!directed && j > i) ann.text.options$col.minus else ann.text.options$col.plus,
							srt=180*(posit.angle(alpha)/pi))
					}
					if (numconn$minus[j,i] > ann.text.options$thresh && directed) {
						txt.num <- if (is.null(numconn.ref)) numconn$minus[j,i]	else numconn$minus[j,i] - numconn.ref$minus[j,i]
						txt.str <- sprintf(paste0("%", if (is.null(numconn.ref)) "" else "+", "1.", ann.text.options$digits, "f"), txt.num)
						text(x=(1-ann.text.options$pos.shift.minus)*x.i.minus+ann.text.options$pos.shift.minus*x.j.minus, 
							y=(1-ann.text.options$pos.shift.minus)*y.i.minus+ann.text.options$pos.shift.minus*y.j.minus, 
							txt.str, 
							cex=ann.text.options$text.cex, 
							col=ann.text.options$col.minus,
							srt=180*(posit.angle(alpha)/pi))
					}
				}
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
					if (numconn$plus[i,i] > ann.text.options$thresh) {
						txt.num <- if (is.null(numconn.ref)) numconn$plus[i,i]	else numconn$plus[i,i] - numconn.ref$plus[i,i]
						txt.str <- sprintf(paste0("%", if (is.null(numconn.ref)) "" else "+", "1.", ann.text.options$digits, "f"), txt.num)
						text(x=cc.plus[round(ann.text.options$pos.shift.plus*nrow(cc.plus)),1], 
							y=cc.plus[round(ann.text.options$pos.shift.plus*nrow(cc.plus)),2], 
							txt.str, 
							cex=ann.text.options$text.cex, 
							col=ann.text.options$col.plus)
					}
					if (numconn$minus[i,i] > ann.text.options$thresh) {
						txt.num <- if (is.null(numconn.ref)) numconn$minus[i,i]	else numconn$minus[i,i] - numconn.ref$minus[i,i]
						txt.str <- sprintf(paste0("%", if (is.null(numconn.ref)) "" else "+", "1.", ann.text.options$digits, "f"), txt.num)						
						text(x=cc.minus[round(ann.text.options$pos.shift.minus*nrow(cc.minus)),1], 
							y=cc.minus[round(ann.text.options$pos.shift.minus*nrow(cc.minus)),2], 
							txt.str, 
							cex=ann.text.options$text.cex, 
							col=ann.text.options$col.minus)
					}
				}
			}
		}
	}
}


plot.Gdiff <- function(mysims, deltaG=NA, xlab="Generation", ylab="Change in G matrix", ylim=c(0,1), ...) {
	gd <- sapply(mysims, function(mysim) {
		out.files  <- list.files(pattern="out.*", path=list.dirs(outdir.all[[mysim]],  full.names=TRUE, recursive=FALSE), full.names=TRUE)
		mean.Gdiff.dyn.cache(out.files, deltaG)
	}, USE.NAMES=TRUE, simplify=FALSE)
	
	gen <-as.numeric(names(gd[[1]]))

	plot(NULL, xlim=c(first.gen, max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	
	for (mysim in mysims) {
		lines(gen, gd[[mysim]], lty=lty.sce[mysim], col=col.sce[mysim])
	}	
}

plot.GPC <- function(mysims, PC=1, xlab="Generation", ylab="Proportion of total variance", ylim=c(0,1), ...) {
	gpc <- sapply(mysims, function(mysim) {
		out.files  <- list.files(pattern="out.*", path=list.dirs(outdir.all[[mysim]],  full.names=TRUE, recursive=FALSE), full.names=TRUE)
		mean.propPC.dyn.cache(out.files, PC)
	}, USE.NAMES=TRUE, simplify=FALSE)
	
	gen <-as.numeric(names(gpc[[1]]))
			
	plot(NULL, xlim=c(first.gen, max(gen)), ylim=ylim, xlab=xlab, ylab=ylab, ...)
	
	for (mysim in mysims) {
		lines(gen, gpc[[mysim]], lty=lty.sce[mysim], col=col.sce[mysim])
	}	
}


plot.Gmat <- function(mysim, gen, absolute=TRUE, cols=NULL) {
	if(is.null(cols))
		cols <- colorRampPalette(col.cor[(if(absolute) 2 else 1):length(col.cor)])(1001)

	st <- if (absolute) 0 else -1
	out.files <- list.files(pattern="out.*", path=list.dirs(outdir.all[[mysim]], full.names=TRUE, recursive=FALSE), full.names=TRUE)
	
	Glist <- Ggen.files.cache(out.files, gen)
	Glist <- lapply(Glist, function(G) if (absolute) abs(cov2cor(G)) else cov2cor(G))
	
	Gmean <- rowMeans(do.call(abind, c(Glist, list(along=3))), dims=2)
	ng <- ncol(Gmean)-1
	
	image(x=1:ng, y=1:ng, t(Gmean[(ng+1):2,2:(ng+1)]), axes=FALSE, xlab="", ylab="", zlim=c(st, 1), col=cols)
}

plot.Gmat.legend <- function(absolute=TRUE, ylab="Genetic correlation", cols=NULL) {
	if(is.null(cols))
		cols <- colorRampPalette(col.cor[(if(absolute) 2 else 1):length(col.cor)])(1001)
	
	st <- if (absolute) 0 else -1
	yl <- paste0(if (absolute) "|" else "", ylab, if (absolute) "|" else "")
	image(y=seq(st, 1, length.out=length(cols)), z=t(as.matrix(seq(st, 1, length.out=length(cols)-1))), zlim=c(st,1), col=cols, xaxt="n", xlab="", ylab=yl)
}