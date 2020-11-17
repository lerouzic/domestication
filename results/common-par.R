# Common parameters (colors, plotting options) shared by all figures.


      
suppressMessages(library(parallel))    # for mclapply
suppressMessages(library(digest))      # to generate cache file names
suppressMessages(library(ellipse))     # plotting ellipses (?)
suppressMessages(library(ade4))        # for mantel.rtest
suppressMessages(library(igraph))      # graph topology calculation
suppressMessages(library(Rcpp))        # fast model
suppressMessages(library(inline))      # fact model

mc.cores <- min(12, detectCores()-1)
cache.dir <- default.cache.dir # defined in ../src/cache.R

window.avg <- 1 # Size of the moving average window
window.norm <- 10 # Size of the window for reaction norms
first.gen  <- 10000  # First generation for the time series
deltaG     <- 1000 # When tracking evolutionary change

connect.threshold <- 0.1
connect.env       <- 0.5

corr.threshold    <- 0.2

col.sel <- c(c="blue", s="blue", p="red", n="black")
lty.sel <- c(c=1, s=1, p=2, n=3)
col.sce <- c(default="black", nobot="black", noselc="black", nosel="gray")
lty.sce <- c(default=1, nobot=2, noselc=3, nosel=3)

col.algo <- c( # modularity algorithms
	walktrap="tomato", 
	fastgreedy="darkolivegreen",
	labelprop="lightblue")
	
col.gl <- c(gain="orange", loss="darkgray")
col.cor <- c("red", "white", "green") # for correlations -1, 0, 1
