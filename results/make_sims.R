#!/usr/bin/env Rscript

# Generate all simulations and run them if not in cache

# A simulation is characterized by a name, a parameter file, and an extended parameter file with bottleneck etc. information

use.cache <- TRUE # Don't run the simulation if already there in the cache

this.dir <- 

param.dir <- file.path(this.dir, "../param")
cache.dir <- file.path(this.dir, "../cache")
src.dir   <- file.path(this.dir, "../src")

source(file.path(src.dir, "makeparam_functions.R"))

all.sims <- rbind(simTest=c("param0-test.txt", "extparam0-test.txt"))

for (sim.name in rownames(all.sim)) {
	
	
	
	
}
