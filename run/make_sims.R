#!/usr/bin/env Rscript

# Generate all simulations and run them if not in cache

# A simulation is characterized by a name, a parameter file, and an extended parameter file with bottleneck etc. information
# Simulations involve the creation of many small files (new parameters every generation due to plasticity). This is largely
# sub-optimal, but difficult to change given the way the simulation program works -- we have to deal with it. 

use.cache <- TRUE # Don't run the simulation if already there in the cache // Not functional
overwrite <- TRUE

options(warn=1)

cl <- commandArgs(trailingOnly = FALSE)
script.dir <- if (grepl (cl[1], pattern="/R")) "." else dirname(strsplit(split="=", cl[grep(cl, pattern="--file")])[[1]][2])

user.dir  <- getwd()

param.dir <- file.path(script.dir, "../param")
cache.dir <- normalizePath(file.path(script.dir, "../cache"))
src.dir   <- file.path(script.dir, "../src")

prog.path <- file.path(src.dir, "Simul_Prog") # This relies on a symbolic link, give the full path otherwise

source(file.path(src.dir, "makeparam_functions.R"))

all.sims <- rbind(
#~ 	simTest    = c("param0-test.txt", "extparam0-test.txt")
	simDefault = c("param0.txt", "extparam0.txt"),
	simNobot   = c("param0.txt", "extparam-nobot.txt"),
	simNoselc  = c("param0.txt", "extparam-noselc.txt"),
	simNosel   = c("param0.txt", "extparam-nosel.txt"),
	simLargenet= c("param-largenet.txt", "extparam-largenet.txt"),
	simSmallsel= c("param-smallsel.txt", "extparam-smallsel.txt"),
	simStrongbot=c("param0.txt", "extparam-strongbot.txt"),
	simNomut   = c("param0.txt", "extparam-nomut.txt"),
	simStrsel  = c("param0.txt", "extparam-strongsel.txt"),
	simIdplast = c("param0.txt", "extparam-idplast.txt"),
	simCstplast= c("param0.txt", "extparam-cstplast.txt")
)

for (sim.name in rownames(all.sims)) {
	cat("Setting up simulation", sim.name, "...\n")
	pars <- create.paramseries(
		file.path(param.dir, all.sims[sim.name, 1]), 
		file.path(param.dir, all.sims[sim.name, 2]), 
		file.path(cache.dir, sim.name), 
		overwrite=overwrite, verbose=TRUE)
	create.launchfile(prog.path, pars$param, pars$out, pars$compressed, file.path(user.dir, paste0(sim.name, "-launch.sh")))
}
