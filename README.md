# Domestication project

Gene network &amp; Domestication project, involving Ewen Burban, Maud Tenaillon, &amp; Arnaud Le Rouzic

This directory contains all the material to run and analyse the simulations.


## Launching simulations

`cd results`
`Rscript ./make_sims.R`

This will generate (i) a  series of parameter files in ../cache, and (ii) a file ./xxx-launch.sh in the current directory. 

It is possible to make a unique launch file for all simulations:

`cat sim* > simAll.sh`

Or to split simulations into batches if too many of them:

`for j in sim*.sh; do for i in {1..10}; do sed -n -e "$(( 100*$i-99 )),$(( 100*$i )) p" -e "$(( 100*$i )) q" $j >> simAll_$i.sh; done; done`

Tar only a subset of the parameters to run on a specific server:

`tar cvfz target.tar.gz t cache/sim*/rep-0{101..200}`

Launching all jobs can be performed with GNU parallel:

`nohup nice parallel --shuf --joblog joblog.log -j 64 -a simAll.sh`

On a slurm server:

```#!/bin/bash
#
#SBATCH -p long                      # partition
#SBATCH --mail-type END
#SBATCH --mail-user toto@domain.com
#SBATCH --cpus-per-task 54

module load parallel tar

parallel --shuf --joblog joblog.log -j $SLURM_CPUS_PER_TASK -a simAll.sh
```

Current simulations:
* simAll_1 -> Taranis
* simAll_2 -> IFBcore
