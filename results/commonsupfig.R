#same than commonfig but for the supplementary simulations
#NB : the simulation list is customisable 

source("../src/analysis_tools.R")

library(parallel)
mc.cores <- min(12, detectCores()-1)

cache.dir <- default.cache.dir # defined in ../src/cache.R

window.avg <- 1 # Size of the moving average window
window.norm <- 10 # Size of the window for reaction norms
first.gen  <- 0  # First generation for the time series
onerep <- function(out.dir) list.dirs(out.dir, full.names=TRUE, recursive=FALSE)[1]

simtag<-c("default","nomut","smallsel","strsel")
simtoplot<-c("simDefault","simNomut","simSmallsel","simStrsel")
tokeep <- "Gen|Mphen|VPhen|FitOpt|MFit|VFit|VarAll"
i<-0
for(I in simtag) {
  i<-i+1;out.dir<-eval(parse(text=paste0("out.dir.",I,"<-file.path(cache.dir,simtoplot[i])")))
  mean.sim.data<-eval(parse(text=paste0("mean.sim.",I,"<-mean.sim.cache(out.dir,colnames.pattern=tokeep)")))
  selpattern<- eval(parse(text=paste0("selpattern.",I,"<-selectionregime.detect(mean.sim.data)")))
  Ndyn<-eval(parse(text=paste0("Ndyn.",I,"<-get.Ndyn.cache(onerep(out.dir))")))
}

col.sel <- c(c="blue", s="blue", p="red", n="black")
lty.sel <- c(c=1, s=1, p=2, n=3)
col.sce <- c(default="black", nobottle="black", noselc="black", nosel="gray")
lty.sce <- c(default=1, nobottle=2, noselc=3, nosel=3)

generation.axis <- function(show.bottleneck=FALSE, ...) {
  mxx <- max(as.numeric(names(Ndyn.default)))
  if (show.bottleneck) {
    bd <- unlist(bottleneck.detect(Ndyn.default))
    toshow <- sort(c(first.gen, mxx, bd))
  } else {
    toshow <- mxx - pretty(c(0, mxx - first.gen), n=3)
    toshow <- toshow[toshow > first.gen]
  }
  axis(1, at=toshow, labels = toshow - mxx, ...)
}

legname <- function(nn) {
  legn <- c(
    default ="Default",
    nobottle="No bottleneck",
    noselc  ="Constant selection")
  
  ifelse(nn %in% names(legn), legn[nn], nn)
}

subpanel <- function(x) {
  title(outer=FALSE,adj=0.025 ,main=x,cex.main=1.4,col="black",line=-1)
}

