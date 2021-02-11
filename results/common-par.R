# Common parameters (colors, plotting options) shared by all figures.
      
suppressMessages(library(parallel))    # for mclapply
suppressMessages(library(digest))      # to generate cache file names
suppressMessages(library(ade4))        # for mantel.rtest
suppressMessages(library(igraph))      # graph topology calculation
suppressMessages(library(Rcpp))        # fast model
suppressMessages(library(inline))      # fact model

source("../src/cache.R")

mc.cores <- min(32, detectCores()-1)
cache.dir <- default.cache.dir # defined in ../src/cache.R

window.avg <- 1 # Size of the moving average window
window.norm <- 10 # Size of the window for reaction norms
first.gen  <- 10000  # First generation for the time series
deltaG     <- 500 # When tracking evolutionary change

connect.threshold <- 0.1
connect.env       <- 0.5

corr.threshold    <- 0.2

cex.legend <- 0.8
panel.width  <- 4
panel.height <- 4
mar.notitle <- c(4.5, 4.5, 0.5, 0.5)

col.sel <- c(c="blue", s="blue", p="red", n="black")
lty.sel <- c(c=1, s=1, p=2, n=3)
col.sce <- c(default="black", nobot="black", noselc="black", nosel="gray")
lty.sce <- c(default=1, nobot=2, noselc=3, nosel=3)

col.algo <- c( # modularity algorithms
	walktrap="tomato", 
	fastgreedy="darkolivegreen",
	labelprop="lightblue")
	
col.gl <- c(Gain="deeppink", Loss="seagreen4")
col.cor <- c("red", "white", "black") # for correlations -1, 0, 1

# Consistent caption text
legname <- function(nn) {
	legn <- c(
		default ="Default",
		nobot   ="No bottleneck",
		noselc  ="No selection switch",
		nosel   ="Drift",
		nomut   ="No new mutations",
		smallsel="Less selected genes",
		strongsel="Strong selection",
		strongbot="Strong bottleneck",
		largenet="Large network", 
		idplast ="Similar plasticity",
		cstplast="Constant plasticity")
		
	ifelse(nn %in% names(legn), legn[nn], nn)
}
