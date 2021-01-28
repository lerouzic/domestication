# Various functions to analyze the network

source("./common-par.R")
source("../src/analysis_tools.R")
source("../src/cache.R")

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


# Extraction of W and G matrices from files

Wgen.files <- function(out.dir, gen) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		if (!gen %in% tt[,"Gen"]) return(NULL)
		W <- tt[tt[,"Gen"] == gen, grepl(colnames(tt), pattern="MeanAll")]
		rm(tt); gc()
		W <- matrix(unlist(W), ncol=sqrt(length(W)), byrow=TRUE)
	}, mc.cores=1)
	ans[!sapply(ans, function(x) length(x) == 1 && is.na(x))]
}

Wgen.files.cache <- function(out.dir, gen) {
	cache.fun(Wgen.files, out.dir=out.dir, gen=gen, cache.subdir="Rcache-Wgen")
}

Wlist.table <- function(out.table) {
	gen <- out.table[,"Gen"]
	out.table <- out.table[, grepl(colnames(out.table), pattern="MeanAll")]
	ans <- lapply(1:nrow(out.table), function(i) {
		W <- out.table[i,]
		matrix(unlist(W), ncol=sqrt(length(W)), byrow=TRUE)
	})
	names(ans) <- as.character(gen)
	ans
}

Wlist.files <- function(out.dir) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		Wlist.table(tt)
		rm(tt); gc()
	}, mc.cores=1)
}

Ggen.files <- function(out.dir, gen, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		if (!gen %in% tt[,"Gen"]) return(NULL)
		cc <- tt[tt[,"Gen"] == gen, grepl(colnames(tt), pattern="CovPhen")]
		rm(tt); gc()
		cc <- matrix(unlist(cc), ncol=sqrt(length(cc)), byrow=TRUE)
		diag(cc)[diag(cc) < 1e-12] <- 1e-12
		cc
	}, mc.cores=1)
	ans[!sapply(ans, function(x) length(x) == 1 && is.na(x))]
}

Ggen.files.cache <- function(out.dir, gen, mc.cores=1) {
	cache.fun(Ggen.files, out.dir=out.dir, gen=gen, mc.cores=mc.cores, cache.subdir="Rcache-Ggen")
}

Glist.table <- function(out.table) {
	gen <- out.table[,"Gen"]
	out.table <- out.table[, grepl(colnames(out.table), pattern="CovPhen")]
	ans <- lapply(1:nrow(out.table), function(i) {
		G <- out.table[i,]
		matrix(unlist(G), ncol=sqrt(length(G)), byrow=TRUE)
	})
	names(ans) <- as.character(gen)
	ans
}

Glist.files <- function(out.dir) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		gl <- Glist.table(tt)
		rm(tt); gc()
		gl
	}, mc.cores=mc.cores)	
}


# Network analysis

cleanW <- function(W, epsilon=connect.threshold, env=connect.env, ...) {
	# ... are additional arguments to modelM2
	if (is.na(epsilon)) epsilon <- sqrt((nrow(W)-1)*0.01^2)
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
	cW <- matrix(0, ncol=ncol(W), nrow=nrow(W))
	cW[distMat > epsilon] <- W[distMat > epsilon]
	cW
}

cleanW.cache <- function(W, epsilon=connect.threshold, env=connect.env, ...) {
	cache.fun(cleanW, W=W, epsilon=epsilon, env=env, ..., cache.subdir="Rcache-cleanW")
}

number.connections <- function(W, ...) {
	cW <- cleanW.cache(W=W, ...)
	sum(cW != 0)
}

number.connections.dyn <- function(out.table) { 
	W.table <- out.table[,grepl(colnames(out.table), pattern="MeanAll")]
	
	net.size <- sqrt(ncol(W.table))
	nb.conn <- sapply(1:nrow(W.table), function(i) { 
			W <- matrix(unlist(W.table[i,]), ncol=net.size, byrow=TRUE)
			if(nrow(W) != ncol(W)) return(NA)
			ans <- try(number.connections(W))
			if (class(ans) == "try-error") NA else ans
		})
	names(nb.conn) <- as.character(out.table[,"Gen"])
	nb.conn
}

inout.connections <- function(W, ...) {
	cW <- cleanW.cache(W=W, ...)
	list(connect.in=rowSums(cW != 0), connect.out=colSums(cW != 0))
}

# in/out connections at a specific generation
inout.gen <- function(out.dir, gen, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		if (!gen %in% tt$Gen) return(NA)
		W <- tt[tt[,"Gen"] == gen, grepl(colnames(tt), pattern="MeanAll")]
		rm(tt); gc()
		W <- matrix(unlist(W), ncol=sqrt(length(W)), byrow=TRUE)
		inout.connections(W)
	}, mc.cores=mc.cores)
	ans[!sapply(ans, function(x) length(x) < 2 || is.na(x))]
}

inout.gen.cache <- function(out.dir, gen, mc.cores=1) {
	cache.fun(inout.gen, out.dir=out.dir, gen=gen, mc.cores=mc.cores, cache.subdir="Rcache-inout")
}

delta.inout <- function(W, Wref) {
	cW <- cleanW.cache(W)
	cWref <- cleanW.cache(Wref)
	
	c(gain=sum(cWref == 0 & cW != 0), loss=sum(cWref != 0 & cW == 0))
}

# Returns two columns: one with the gains, 
delta.inout.dyn <- function(out.table, deltaG=1, mc.cores=1) {
	gen <- out.table[,"Gen"]
	seqgen <- seq(1, length(gen), length.out=min(length(gen), 1+gen[length(gen)] %/% deltaG))
	listW <- Wlist.table(out.table[seqgen,])
	
	ans <- mclapply(1:(length(listW)-1), function(i) delta.inout(listW[[i+1]], listW[[i]]), mc.cores=mc.cores)
	ans <- do.call(rbind, ans)
	rownames(ans) <- as.character(gen[seqgen])[-1]
	ans
}

mean.delta.inout.dyn <- function(out.dir, deltaG=NA, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		delta.inout.dyn(tt, deltaG, mc.cores=1)
	}, mc.cores=mc.cores)
	l.ans <- sapply(ans, nrow)
	ans <- ans[l.ans == max(l.ans)] # removing incomplete data sets (ongoing simulations)
	aa <- do.call(abind, c(ans, list(along=3)))
	rowMeans(aa, dims=2)
}

mean.delta.inout.dyn.cache <- function(out.dir, deltaG=NA, mc.cores=1) {
	cache.fun(mean.delta.inout.dyn, out.dir=out.dir, deltaG=deltaG, mc.cores=mc.cores, cache.subdir="Rcache-dinout")
}


delta.Wdiff <- function(W, W.ref) {
	# Euclidian distance line by line
	sqrt(rowSums((W-W.ref)^2))
}

delta.Wdiff.dyn <- function(out.table, deltaG=NA, mc.cores=1) {
	gen <- out.table[,"Gen"]
	seqgen <- seq(1, length(gen), length.out=min(length(gen), 1+gen[length(gen)] %/% deltaG))
	listW <- Wlist.table(out.table[seqgen,])
	
	ans <- mclapply(1:(length(listW)-1), function(i) delta.Wdiff(listW[[i+1]], listW[[i]]), mc.cores=mc.cores)
	names(ans) <- as.character(gen[seqgen])[-1]
	do.call(rbind, ans)		
}

mean.Wdiff.dyn <- function(out.dir, deltaG=NA, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		delta.Wdiff.dyn(tt, deltaG, mc.cores=1)
	}, mc.cores=mc.cores)
	
	arr <- do.call(abind, c(ans, list(along=3)))
	rowMeans(arr, dims=2)
}

mean.Wdiff.dyn.cache <- function(out.dir, deltaG=NA, mc.cores=1) {
	cache.fun(mean.Wdiff.dyn, out.dir=out.dir, deltaG=deltaG, mc.cores=mc.cores, cache.subdir="Rcache-Wdiff")
}


delta.Gdiff <- function(G, G.ref) {
	# A bit of cleaning is necessary : first generation and "environmental" gene can mess up the vcov
	if (any(!is.finite(G)) || any(!is.finite(G.ref))) return(NA)
	diag(G)[diag(G) < 1e-8] <- 1e-8
	diag(G.ref)[diag(G.ref) < 1e-8] <- 1e-8
	# covtransf is defined in common-precalc.R, it turns vcov into distance matrices
	mt <- try(mantel.rtest(covtransf(G), covtransf(G.ref), nrepet=1)$obs)
	if (class(mt) == "try-error") mt <- 1
	1 - mt
}

delta.Gdiff.dyn <- function(out.table, deltaG=NA, mc.cores=1) {
	gen <- out.table[,"Gen"]
	seqgen <- seq(1, length(gen), length.out=min(length(gen), 1+gen[length(gen)] %/% deltaG))
	listG <- Glist.table(out.table[seqgen,])
	
	ans <- mclapply(1:(length(listG)-1), function(i) delta.Gdiff(listG[[i+1]], listG[[i]]), mc.cores=mc.cores)
	ans <- do.call(c, ans)
	names(ans) <- as.character(gen[seqgen])[-1]
	ans	
}

mean.Gdiff.dyn <- function(out.dir, deltaG=NA, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		delta.Gdiff.dyn(tt, deltaG, mc.cores=1)
	}, mc.cores=mc.cores)

	ansl <- sapply(ans, length)
	ans <- ans[ansl == max(ansl)]
	aa <- do.call(rbind, ans)
	colMeans(aa, na.rm=TRUE)
}

mean.Gdiff.dyn.cache <- function(out.dir, deltaG=NA, mc.cores=1) {
	cache.fun(mean.Gdiff.dyn, out.dir=out.dir, deltaG=deltaG, mc.cores=mc.cores, cache.subdir="Rcache-Gdiff")
}

Gcor.dyn <- function(out.table) {
	.meancor <- function(G) {
		diag(G)[diag(G) < 1e-12] <- 1e-12 # Prevent numerical errors 
		C <- cov2cor(G)
		diag(C) <- NA # We don't want to count the diagonal when computing the mean
		mean(abs(C), na.rm=TRUE)
	}
	Glist <- Glist.table(out.table)
	sapply(Glist, .meancor)
}


mean.Gcor.dyn.files <- function(out.dir, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		gc <- Gcor.dyn(tt)
		rm(tt)
		gc
	}, mc.cores=mc.cores)
	
	ansl <- sapply(ans, length)
	ans <- ans[ansl==max(ansl)]
	aa <- do.call(rbind, ans)
	
	colMeans(aa, na.rm=TRUE)
}

mean.Gcor.dyn.files.cache <- function(out.dir, mc.cores=1) {
	cache.fun(mean.Gcor.dyn.files, out.dir=out.dir, mc.cores=mc.cores, cache.subdir="Rcache-Gcor")
}

WFUN.dyn <- function(out.table, WFUN="mean", deltaG=1) {
	gen <- out.table[,"Gen"]
	seqgen <- seq(1, length(gen), length.out=min(length(gen), 1+gen[length(gen)] %/% deltaG))
	listW <- Wlist.table(out.table[seqgen,])
	ww <- unlist(lapply(listW, eval(parse(text=WFUN))))
	ww
}

WFUN.dyn.files <- function(out.dir, WFUN="mean", deltaG=1, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		wm <- WFUN.dyn(tt, WFUN=WFUN, deltaG=deltaG)
		rm(tt)
		wm
	}, mc.cores=mc.cores)
	ansl <- sapply(ans, length)
	ans <- ans[ansl==max(ansl)]
	aa <- do.call(rbind, ans)
	colMeans(aa, na.rm=TRUE)
}

WFUN.dyn.files.cache <- function(out.dir, WFUN="mean", deltaG=1, mc.cores=1) {
	cache.fun(WFUN.dyn.files, out.dir=out.dir, WFUN=WFUN, deltaG=deltaG, mc.cores=mc.cores, cache.subdir=paste0("Rcache-W", WFUN))
}


#~ netf.dyn <- function(out.table, what, directed=TRUE, deltaG=NA) {
#~ 	gen <- out.table[,"Gen"]
#~ 	seqgen <- seq(1, length(gen), length.out=min(length(gen), 1+gen[length(gen)] %/% deltaG))
#~ 	listW <- Wlist.table(out.table[seqgen,])
#~ 	ww <- unlist(lapply(listW, eval(parse(text=WFUN))))
#~ 	ww
#~ }

#~ netf.dyn.files <- function(outdir, what, directed=TRUE, deltaG=NA, mc.cores=1) {
#~ 	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
#~ 	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
#~ 	ans <- mclapply(out.files, function(ff) {
#~ 		tt <- read.table(ff, header=TRUE)
#~ 		wm <- netf.dyn(tt, what=what, directed=directed, deltaG=deltaG)
#~ 		rm(tt)
#~ 		wm
#~ 	}, mc.cores=mc.cores)
#~ 	ansl <- sapply(ans, length)
#~ 	ans <- ans[ansl==max(ansl)]
#~ 	aa <- do.call(rbind, ans)
#~ 	colMeans(aa, na.rm=TRUE)
#~ }

#~ netf.dyn.files.cache <- function(outdir, what, directed=TRUE, deltaG=NA, mc.cores=1) {
#~ 	cache.fun(netf.dyn.files, out.dir=out.dir, what=what, directed=directed, deltaG=deltaG, mc.cores=mc.cores, cache.subdir="Rcache-netwf")
#~ }



propPC.dyn <- function(out.table, PC=1, mc.cores=1) {
	gen <- out.table[,"Gen"]
	listG <- Glist.table(out.table)
	
	ans <- mclapply(listG, function(G) {
		ee <- eigen(G)
		ee$values[PC]/sum(ee$values)
	}, mc.cores=mc.cores)
	aa <- do.call(c, ans)
	names(aa) <- as.character(gen)
	aa
}

mean.propPC.dyn <- function(out.files, PC=1, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		propPC.dyn(tt, PC=PC, mc.cores=1)
	}, mc.cores=mc.cores)
	ansl <- sapply(ans, length)
	ans <- ans[ansl == max(ansl)]
	aa <- do.call(rbind, ans)
	colMeans(aa)
}

mean.propPC.dyn.cache <- function(out.dir, PC=1, mc.cores=1) {
	cache.fun(mean.propPC.dyn, out.dir=out.dir, PC=PC, mc.cores=mc.cores, cache.subdir="Rcache-propPC")
}

erankG.dyn <- function(out.table, mc.cores=1) {
	erank <- function(M) {
		# calculation of the effective rank as defined in 
		# Roy, O., & Vetterli, M. (2007). 
		# The effective rank: A measure of effective dimensionality. 
		# In 2007 15th European Signal Processing Conference (pp. 606-610). IEEE.
		ee <- abs(eigen(M)$values)
		p  <- ee/sum(ee)
		sa <- - sum(ifelse(p==0, 0, p*log(p))) # Shannon entropy
		exp(sa)
	}
	
	gen <- out.table[,"Gen"]
	listG <- Glist.table(out.table)
	
	ans <- mclapply(listG, function(G) {
		erank(G)
	}, mc.cores=mc.cores)
	aa <- do.call(c, ans)
	names(aa) <- as.character(gen)
	aa
}


mean.erankG.dyn <- function(out.dir, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	ans <- mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		erankG.dyn(tt, mc.cores=1)
	}, mc.cores=mc.cores)
	ansl <- sapply(ans, length)
	ans <- ans[ansl == max(ansl)]
	aa <- do.call(rbind, ans)
	colMeans(aa)
}

mean.erankG.dyn.cache <- function(out.dir, mc.cores=1) {
	cache.fun(mean.erankG.dyn, out.dir=out.dir, mc.cores=mc.cores, cache.subdir="Rcache-erankG")
}


# Average out all network connections from a directory 
mean.connect <- function(out.dir, max.reps=Inf, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	tt <- results.table(out.files, mc.cores, max.reps)
	nn <- mclapply(tt, number.connections.dyn, mc.cores=mc.cores)
	ans <- colMeans(do.call(rbind, nn), na.rm=TRUE)
	rm(tt)
	gc()
	return(ans)
}

mean.connect.cache <- function(out.dir, max.reps=Inf, mc.cores=1) {
	cache.fun(mean.connect, out.dir=out.dir, max.reps=max.reps, mc.cores=mc.cores, cache.subdir="Rcache-connect")
}

# Returns a list of complex community objects according to several igraph algorithms
communities <- function(W, directed=FALSE, ...) {
	cW <- cleanW.cache(W, ...)
	
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
communities.dyn <- function(out.table, directed=FALSE, mc.cores=1) {
	W.table <- out.table[,grepl(colnames(out.table), pattern="MeanAll")]
	
	net.size <- sqrt(ncol(W.table))
	comm <- mclapply(1:nrow(W.table), function(i) { 
			W <- matrix(unlist(W.table[i,]), ncol=net.size, byrow=TRUE)
			stopifnot(nrow(W) == ncol(W))
			communities(W, directed=directed)
		}, mc.cores=mc.cores)
	names(comm) <- as.character(out.table[,"Gen"])
	comm
}

communities.dyn.files <- function(out.dir, directed=FALSE, mc.cores=1) {
	out.reps <- list.dirs(out.dir, full.names=TRUE, recursive=FALSE)
	out.files <- list.files(pattern="out.*", path=out.reps, full.names=TRUE)
	 mclapply(out.files, function(ff) {
		tt <- read.table(ff, header=TRUE)
		cc <- communities.dyn(tt, directed=directed, mc.cores=1)
		rm(tt); gc()
		return(cc)
		}, mc.cores=mc.cores)
}

communities.dyn.files.cache <- function(out.dir, directed=FALSE, mc.cores=1) {
	cache.fun(communities.dyn.files, out.dir=out.dir, directed=directed, mc.cores=mc.cores, cache.subdir="Rcache-commdyn")
}

numconn.groups <- function(W, groups, ...) {
	cW <- cleanW.cache(W=W, ...)
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

mean.numconn.groups <- function(listW, groups, count.diag=NA, mc.cores=1, ...) {
	if (is.na(count.diag)) count.diag <- sum(sapply(listW, function(W) sum(diag(W)!=0))) != 0
	all.nconn <- mclapply(listW, function(W) numconn.groups(W=W, groups=groups, ...), mc.cores=mc.cores)
	tg <- table(groups)
	norm <- tg %*% t(tg)
	if (!count.diag)
		diag(norm) <- diag(norm)-tg
	tplus <- do.call(abind, c(lapply(all.nconn, function(x) x$plus), list(along=3)))
	tminus <-  do.call(abind, c(lapply(all.nconn, function(x) x$minus), list(along=3)))
	list(plus=rowMeans(tplus, dims=2)/norm, minus=rowMeans(tminus, dims=2)/norm)
}


numcorrgen.groups <- function(R, groups, cutoff) {
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


mean.numcorrgen.groups <- function(listR, groups, cutoff=corr.threshold, remove.e=TRUE, mc.cores=1) {
	all.ncorr <- mclapply(listR, function(R) numcorrgen.groups(R=R, groups=groups, cutoff=cutoff), mc.cores=mc.cores)
	tg <- table(groups)
	norm <- tg %*% t(tg)
	diag(norm) <- diag(norm)-tg
	if (remove.e) norm <- norm[rownames(norm)!="e",colnames(norm)!="e"]
	tplus <- do.call(abind, c(lapply(all.ncorr, function(x) if (remove.e) x$plus[rownames(x$plus)!="e",colnames(x$plus)!="e"] else x$plus), list(along=3)))
	tminus <-  do.call(abind, c(lapply(all.ncorr, function(x) if (remove.e) x$minus[rownames(x$minus)!="e",colnames(x$minus)!="e"] else x$minus), list(along=3)))
	list(plus=rowMeans(tplus, dims=2)/norm, minus=rowMeans(tminus, dims=2)/norm)
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

mean.numcorrenv.groups <- function(listW, groups, cutoff=0.1, mc.cores=1) {
	all.ncorr <- mclapply(listW, function(W) numcorrenv.groups(W=W, groups=groups, cutoff=cutoff), mc.cores=mc.cores)
	tg <- table(groups)
	norm <- tg %*% t(tg)
	diag(norm) <- diag(norm)-tg
	tplus <- do.call(abind, c(lapply(all.ncorr, function(x) x$plus), list(along=3)))
	tminus <-  do.call(abind, c(lapply(all.ncorr, function(x) x$minus), list(along=3)))
	list(plus=rowMeans(tplus, dims=2)/norm, minus=rowMeans(tminus, dims=2)/norm)
}
