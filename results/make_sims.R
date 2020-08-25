#!/usr/bin/env Rscript

# Generate all simulations and run them if not in cache

# A simulation is characterized by a name, a parameter file, and an extended parameter file with bottleneck etc. information

use.cache <- TRUE # Don't run the simulation if already there in the cache
overwrite <- TRUE

cl <- commandArgs(trailingOnly = FALSE)
script.dir <- dirname(strsplit(split="=", cl[grep(cl, pattern="--file")])[[1]][2])
user.dir  <- getwd()

param.dir <- file.path(script.dir, "../param")
cache.dir <- file.path(script.dir, "../cache")
src.dir   <- file.path(script.dir, "../src")

prog.path <- file.path(src.dir, "Simul_Prog") # This relies on a symbolic link, give the full path otherwise

source(file.path(src.dir, "makeparam_functions.R"))

all.sims <- rbind(simTest=c("param0-test.txt", "extparam0-test.txt"))

for (sim.name in rownames(all.sims)) {
	pars <- create.paramseries(
		file.path(param.dir, all.sims[sim.name, 1]), 
		file.path(param.dir, all.sims[sim.name, 2]), 
		file.path(cache.dir, sim.name), 
		overwrite=overwrite)
	create.launchfile(prog.path, pars$param, pars$out, file.path(user.dir, paste0(sim.name, "-launch.sh")))
}
