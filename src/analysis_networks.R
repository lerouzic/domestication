# Various functions to analyze the network


source("../src/analysis_tools.R")
source("../src/cache.R")

suppressMessages(library(igraph))
suppressMessages(library(Rcpp))
suppressMessages(library(inline))
suppressMessages(library(digest))
suppressMessages(library(ellipse))

cppFunction('
	List internal_loop_cpp(const NumericMatrix &W, const NumericVector &S0, double a, double env, unsigned int steps, unsigned int measure) {
		double lambda = (1-a)/a;
		double mu     = 1/(a*(1-a));
		NumericMatrix sto (S0.size(), steps+1);
		NumericVector sumx (S0.size());
		NumericVector sumx2 (S0.size()); 
		sto(0,0) = env;
		for (unsigned int i = 1; i < S0.size(); i++)
			sto(i,0) = S0(i);
		for (unsigned int t = 1; t <= steps; t++) {
			for (unsigned int i = 0; i < S0.size(); i++) {
				double tmp = 0.;
				if (i == 0) { // The signal gene
					tmp = env;
				} else {
					for (unsigned int j = 0; j < S0.size(); j++) {
						tmp += sto(j,t-1) * W(i,j);
					}
					tmp =  1. / (1. + lambda * exp(-mu*tmp));
				}
				
				sto(i,t) = tmp;
				
				if (t > steps-measure) {
					sumx(i) += tmp;
					sumx2(i) += tmp*tmp;
				}
			}
		}
		for (unsigned int i = 0; i < S0.size(); i++) {
			sumx(i) /= static_cast<double>(measure); // sumx(i) now contains the mean
			sumx2(i) /= static_cast<double>(measure);
			sumx2(i) += -sumx(i)*sumx(i); // sumx2(i) now contains the variance
		}
		return List::create(Named("full")=sto, Named("mean")=sumx, Named("var")=sumx2);
	}')


model.M2 <- function(W, S0=rep(a, nrow(W)), a=0.2, env=0.5, steps=20, measure=4, full=FALSE) {
    ans <- internal_loop_cpp(W, S0, a, env, steps, measure)
	if (!full) ans$full <- NULL
	return(ans)
}

cleanW <- function(W, epsilon=NULL, env=0.5, ...) {
	# ... are additional arguments to modelM2

	if (is.null(epsilon)) epsilon <- sqrt((nrow(W)-1)*0.01^2)
	distMat <- matrix(0, ncol=ncol(W), nrow=nrow(W))
	ref.expr <- model.M2(W, env=env, ...)$mean
	
	for (i in 1:nrow(W))
		for (j in 1:ncol(W)) {
			myW <- W
			myW[i,j] <- 0
			new.expr <- model.M2(myW, env=env, ...)$mean
			dd <- sqrt(sum((ref.expr[-1]-new.expr[-1])^2))
			distMat[i,j] <- dd
		}
	cleanW <- matrix(0, ncol=ncol(W), nrow=nrow(W))
	cleanW[distMat > epsilon] <- W[distMat > epsilon]
	cleanW
}

cleanW.cache <- function(W, epsilon=NULL, env=0.5, ...) {
	cache.fun(cleanW, W=W, epsilon=epsilon, env=env, ..., cache.subdir="cleanW")
}

number.connections <- function(W, epsilon=NULL, env=0.5, ...) {
	cW <- cleanW.cache(W=W, epsilon=epsilon, env=env, ...)
	sum(cW != 0)
}

number.connections.dyn <- function(out.table, epsilon=NULL, env=0.5) { # if env == NULL, MPhen1 is used instead
	W.table <- out.table[,grepl(colnames(out.table), pattern="MeanAll")]
	env <- if(is.null(env)) out.table[,"MPhen1"] else rep(env, nrow(out.table))
	
	net.size <- sqrt(ncol(W.table))
	nb.conn <- sapply(1:nrow(W.table), function(i) { 
			W <- matrix(unlist(W.table[i,]), ncol=net.size, byrow=TRUE)
			stopifnot(nrow(W) == ncol(W))
			number.connections(W, epsilon=epsilon, env=env[i])
		})
	names(nb.conn) <- as.character(out.table[,"Gen"])
	nb.conn
}

inout.connections <- function(W, epsilon=NULL, env=0.5, ...) {
	cW <- cleanW.cache(W=W, epsilon=epsilon, env=env, ...)
	list(connect.in=rowSums(cW != 0), connect.out=colSums(cW != 0))
}

# Average out all network connections from a directory 
mean.connect <- function(out.dir, env=0.5, epsilon=NULL, max.reps=Inf, mc.cores=detectCores()-1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	tt <- results.table(out.files, mc.cores, max.reps)
	nn <- mclapply(tt, number.connections.dyn, env=env, epsilon=epsilon, mc.cores=mc.cores)
	ans <- colMeans(do.call(rbind, nn))
	rm(tt)
	gc()
	return(ans)
}

mean.connect.cache <- function(out.dir, env=0.5, epsilon=NULL, max.reps=Inf, mc.cores=detectCores()-1) {
	cache.fun(mean.connect, out.dir=out.dir, env=env, epsilon=epsilon, max.reps=max.reps, mc.cores=detectCores()-1, cache.subdir="connect")
}

# Returns a list of complex community objects according to several igraph algorithms
communities <- function(W, epsilon=NULL, env=0.5, directed=FALSE, ...) {
	cW <- cleanW.cache(W, epsilon=epsilon, env=env, ...)
	
	Wgraph <- igraph::graph_from_adjacency_matrix(sign(cW))
	if (!directed) 
		Wgraph <- igraph::as.undirected(Wgraph)

	list(
#~ 		edge     = igraph::edge.betweenness.community(Wgraph),
		walktrap = igraph::walktrap.community(Wgraph),
		fastgreedy= igraph::fastgreedy.community(Wgraph),
		labelprop= igraph::label.propagation.community(Wgraph))
}

# Returns the communities algorithm for all generations of a data dable. Beware, the returned object is complex and needs to be further processed
communities.dyn <- function(out.table, epsilon=NULL, env=0.5, directed=FALSE, mc.cores=1) { # if env == NULL, MPhen1 is used instead
	W.table <- out.table[,grepl(colnames(out.table), pattern="MeanAll")]
	env <- if(is.null(env)) out.table[,"MPhen1"] else rep(env, nrow(out.table))
	
	net.size <- sqrt(ncol(W.table))
	comm <- mclapply(1:nrow(W.table), function(i) { 
			W <- matrix(unlist(W.table[i,]), ncol=net.size, byrow=TRUE)
			stopifnot(nrow(W) == ncol(W))
			communities(W, epsilon=epsilon, env=env[i], directed=directed)
		}, mc.cores=mc.cores)
	names(comm) <- as.character(out.table[,"Gen"])
	comm
}

communities.dyn.cache <- function(out.table, epsilon=NULL, env=0.5, directed=FALSE, mc.cores=1) {
	cache.fun(communities.dyn, out.table=out.table, epsilon=epsilon, env=env, directed=directed, mc.cores=mc.cores)
}

numconn.groups <- function(W, groups, epsilon=NULL, env=0.5, ...) {
	cW <- cleanW.cache(W=W, epsilon=epsilon, env=env, ...)
	ug <- sort(unique(groups))
	nconn.plus <- nconn.minus <- matrix(0, ncol=length(ug), nrow=length(ug))
	rownames(nconn.plus) <- rownames(nconn.minus) <- colnames(nconn.plus) <- colnames(nconn.minus) <- ug
	# Very slow double for loop
	for (i in 1:nrow(cW))
		for (j in 1:ncol(cW)) {
			if (cW[i,j] > 0) nconn.plus[groups[i],groups[j]] <- nconn.plus[groups[i],groups[j]]+1
			if (cW[i,j] < 0) nconn.minus[groups[i],groups[j]] <- nconn.minus[groups[i],groups[j]]+1
		}
	list(plus=nconn.plus, minus=nconn.minus)
}

mean.numconn.groups <- function(listW, groups, epsilon=NULL, env=0.5, count.diag=NA, mc.cores=detectCores()-1, ...) {
	if (is.na(count.diag)) count.diag <- sum(sapply(listW, function(W) sum(diag(W)!=0))) != 0
	all.nconn <- mclapply(listW, function(W) numconn.groups(W=W, groups=groups, epsilon=epsilon, env=env, ...), mc.cores=mc.cores)
	tg <- table(groups)
	norm <- tg %*% t(tg)
	if (!count.diag)
		diag(norm) <- diag(norm)-tg
	tplus <- do.call(abind, c(lapply(all.nconn, function(x) x$plus), list(along=3)))
	tminus <-  do.call(abind, c(lapply(all.nconn, function(x) x$minus), list(along=3)))
	list(plus=rowMeans(tplus, dims=2)/norm, minus=rowMeans(tminus, dims=2)/norm)
}


numcorrgen.groups <- function(R, groups, cutoff=0.1) {
	cR <- R
	cR[abs(cR) < 1e-10] <- 0
	diag(cR)[diag(cR) == 0] <- 1 # It happens with the environment
	cR <- cov2cor(cR)
	cR[abs(cR) < cutoff] <- 0
	diag(cR) <- 0
	ug <- sort(unique(groups))
	nconn.plus <- nconn.minus <- matrix(0, ncol=length(ug), nrow=length(ug))
	rownames(nconn.plus) <- rownames(nconn.minus) <- colnames(nconn.plus) <- colnames(nconn.minus) <- ug
	# Very slow double for loop
	for (i in 1:nrow(cR))
		for (j in 1:ncol(cR)) {
			if (cR[i,j] > 0) nconn.plus[groups[i],groups[j]] <- nconn.plus[groups[i],groups[j]]+1
			if (cR[i,j] < 0) nconn.minus[groups[i],groups[j]] <- nconn.minus[groups[i],groups[j]]+1
		}
	list(plus=nconn.plus, minus=nconn.minus)
}


mean.numcorrgen.groups <- function(listR, groups, cutoff=0.1, remove.e=TRUE, mc.cores=detectCores()-1) {
	all.ncorr <- mclapply(listR, function(R) numcorrgen.groups(R=R, groups=groups, cutoff=cutoff), mc.cores=mc.cores)
	tg <- table(groups)
	norm <- tg %*% t(tg)
	diag(norm) <- diag(norm)-tg
	if (remove.e) norm <- norm[rownames(norm)!="e",colnames(norm)!="e"]
	tplus <- do.call(abind, c(lapply(all.ncorr, function(x) if (remove.e) x$plus[rownames(x$plus)!="e",colnames(x$plus)!="e"] else x$plus), list(along=3)))
	tminus <-  do.call(abind, c(lapply(all.ncorr, function(x) if (remove.e) x$minus[rownames(x$minus)!="e",colnames(x$minus)!="e"] else x$minus), list(along=3)))
	list(plus=rowMeans(tplus, dims=2)/norm, minus=rowMeans(tminus, dims=2)/norm)
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


# Extraction of G matrices from files

Ggen <- function(files, gen) {
	ans <- mclapply(files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		if (!gen %in% tt[,"Gen"]) return(NA)
		cc <- tt[tt[,"Gen"] == gen, grepl(colnames(tt), pattern="CovPhen")]
		rm(tt); gc()
		cc <- matrix(unlist(cc), ncol=sqrt(length(cc)), byrow=TRUE)
		diag(cc)[diag(cc) < 1e-6] <- 1
		cc
	}, mc.cores=1)
	ans[!sapply(ans, function(x) length(x) == 1 && is.na(x))]
}

Ggen.cache <- function(files, gen) {
	cache.fun(Ggen, files=files, gen=gen, cache.subdir="Ggen")
}



##### Tools for coexpression networks (not used -> deletion?)

scalefree.regression <- function(M) {
	# Returns the (adjusted) r^2 of the scale free regresssion of the number of connections
	stopifnot(is.matrix(M), ncol(M) == nrow(M), isSymmetric(M))
	diag(M) <- 0
	nconn <- apply(M, 1, function(x) sum(x != 0))
	freqtt <- table(nconn)/length(nconn)
	freqtt <- freqtt[names(freqtt) != "0"]
	ll <- lm(log(freqtt) ~ log(as.numeric(names(freqtt))))
	summary(ll)$adj.r.squared
}

cor.pval <- function(R, N, p.adjust.method="holm") {
	stopifnot(is.matrix(R), nrow(R) == ncol(R), isSymmetric(R), all(R >= -1), all(R <= 1), N > 2)
	t.matrix <- abs(R) * sqrt(N-2)/sqrt(1-R*R)
	p.matrix <- 2*(1-pt(t.matrix, N-2))
	p.corr <- p.adjust(p.matrix[upper.tri(p.matrix)], method=p.adjust.method)
	# Now ensure the output is symmetric
	p.matrix[upper.tri(p.matrix)] <- p.corr
	p.matrix <- t(p.matrix)
	p.matrix[upper.tri(p.matrix)] <- p.corr
	p.matrix
}


numcorrenv.groups <- function(W, groups, cutoff=0.1, envs = seq(0.01, 0.99, length.out=21), ...) {
	dd <- do.call(rbind, lapply(envs, function(env) model.M2(W=W, env=env, ...)$mean))
	cR <- var(dd)
	cR[abs(cR) < 1e-10] <- 0
	
	cR <- cR/diag(cR)
	cR[abs(cR) < cutoff] <- 0
	diag(cR) <- 0
	ug <- sort(unique(groups))
	nconn.plus <- nconn.minus <- matrix(0, ncol=length(ug), nrow=length(ug))
	rownames(nconn.plus) <- rownames(nconn.minus) <- colnames(nconn.plus) <- colnames(nconn.minus) <- ug
	# Very slow double for loop
	for (i in 1:nrow(cR))
		for (j in 1:ncol(cR)) {
			if (cR[i,j] > 0) nconn.plus[groups[i],groups[j]] <- nconn.plus[groups[i],groups[j]]+1
			if (cR[i,j] < 0) nconn.minus[groups[i],groups[j]] <- nconn.minus[groups[i],groups[j]]+1
		}
	list(plus=nconn.plus, minus=nconn.minus)
}

mean.numcorrenv.groups <- function(listW, groups, cutoff=0.1, mc.cores=detectCores()-1) {
	all.ncorr <- mclapply(listW, function(W) numcorrenv.groups(W=W, groups=groups, cutoff=cutoff), mc.cores=mc.cores)
	tg <- table(groups)
	norm <- tg %*% t(tg)
	diag(norm) <- diag(norm)-tg
	tplus <- do.call(abind, c(lapply(all.ncorr, function(x) x$plus), list(along=3)))
	tminus <-  do.call(abind, c(lapply(all.ncorr, function(x) x$minus), list(along=3)))
	list(plus=rowMeans(tplus, dims=2)/norm, minus=rowMeans(tminus, dims=2)/norm)
}
