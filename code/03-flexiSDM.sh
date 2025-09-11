#!/bin/bash

#SBATCH --job-name=03flexiSDM
#SBATCH --nodes=1
#SBATCH --cpus-per-task=12
#SBATCH --mem-per-cpu=40000
#SBATCH --hint=nomultithread
#SBATCH --partition=cpu
#SBATCH --account=eesc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cscher@usgs.gov
#SBATCH -o %x-%j.out
#SBATCH -e %x-%j.err

#mem-per-cpu was 30000

## Load modules (you can see options using 'module avail')
module load cray-R-spatial


## Set the name of the R script to run, and the directory in which to save outputs
script=/caldera/hovenweep/projects/usgs/ecosystems/eesc/cscher/iSDM-framework/code/03-flexiSDM.R


# Inputs
num=$1
block=$2
maxchain=$3
local=$4

# run your script with inputs
srun --cpu-bind=none Rscript $script $num $block $maxchain $local
