# Domestication project

Gene network &amp; Domestication project, involving Ewen Burban, Maud Tenaillon, &amp; Arnaud Le Rouzic

This directory contains all the material to run and analyse the simulations.

Structure of the sub-directories:
* **param**: parameter files for the simulations
* **src**  : R helper source files for various calculations + symbolic link to the C++ simulation software
* **run**  : Rscripts for running simulations
* **results**: R scripts for figures (require simulation results)
* **cache**: cache directory, contains the simulation results (cache/sim*) + results of heavy calculations when generating figures (cache/Rcache*)

Dependencies 
* The simulation program Simul_Prog and its dependencies
* A recent version of R (> 3.4)
* The following R libraries:
** parallel
** digest
** ade4
** igraph
** Rcpp
** inline

## Launching simulations

`cd run`
`Rscript ./make_sims.R`

This will generate (i) a  series of parameter files in ../cache, and (ii) a file ./xxx-launch.sh in the current directory. 

Launching all jobs can be performed with GNU parallel (or any other task manager):

`nohup nice parallel --shuf --joblog joblog.log -j 64 -a simAll.sh`

## Generating figures

Each figure has an associated R script : figX.R, which generates figX.pdf. Run them individually (`Rscript figX.R`) or all at once (`for i in fig*.R; do Rscript $i; done`. Running figures the first time should be slow, and much faster the next time because of the storage of intermediate results in the cache directory. 
