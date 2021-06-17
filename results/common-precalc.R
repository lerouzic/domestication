
source("./common-par.R")

source("../src/analysis_tools.R")
source("../src/analysis_networks.R")

cor2dist <- function(r) as.dist(sqrt(2*abs(1-r)))
covtransf <- function(x) cor2dist(cov2cor(x))


onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]
tokeep <- "Gen|Mphen|VPhen|FitOpt|MFit|VFit|VarAll"

simtag <- c(
	default="simDefault", 
	nobot  ="simNobot", 
	noselc ="simNoselc",
	nosel  ="simNosel",
	nomut  ="simNomut",
	smallsel="simSmallsel",
	strongsel="simStrsel",
	strongbot="simStrongbot",
	largenet="simLargenet", 
	idplast ="simIdplast", 
	cstplast="simCstplast",
	africe  ="simAfrice",
	pemil   ="simPemil",
	tomato  ="simTomato",
	step8   ="simStep8",
	step24  ="simStep24",
	stab2   ="simStab2",
	nostab  ="simNostab"
	)
	
outdir.all     <- lapply(simtag, function(nn) file.path(cache.dir, nn))
meansim.all    <- lapply(outdir.all, function(od) mean.sim.cache(od, colnames.pattern=tokeep, mc.cores=mc.cores))
selpattern.all <- lapply(meansim.all, function(ms) selectionregime.detect(ms)[-1])
Ndyn.all       <- lapply(outdir.all, function(od) get.Ndyn.cache(onerep(od)))

