#!/usr/bin/env Rscript

# Figure B: evolution of molecular variation during domestication

# Four panels: Default, no bottleneck, no change in selection, no selection. 

library(parallel)
mc.cores <- min(12, detectCores()-1)

source("../src/analysis_tools.R")

out.dir.default <- "../cache/simDefault"
out.dir.nobottle <- "../cache/simNobot"
out.dir.noselc   <- "../cache/simNoselc"
out.dir.nosel    <- "../cache/simNosel"

