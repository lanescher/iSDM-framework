#!/bin/bash

#SBATCH --job-name=01flexiSDM
#SBATCH --array=1-4:1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem-per-cpu=10000
#SBATCH --hint=nomultithread
#SBATCH --time=00-6:00:00
#SBATCH --partition=cpu
#SBATCH --account=eesc
#SBATCH --mail-type=ALL
#SBATCH --mail-user=rmummah@usgs.gov
#SBATCH -o %x-%A-%a.out
#SBATCH -e %x-%A-%a.err

#mem-per-cpu was 30000

## Load modules (you can see options using 'module avail')
module load cray-R-spatial


## Set the name of the R script to run, and the directory in which to save outputs
script=/caldera/hovenweep/projects/usgs/ecosystems/eesc/rmummah/iSDM-framework/code/01-flexiSDM.R


# Inputs
num=$1
# block = array number $SLURM_ARRAY_TASK_ID
local=$2

# run your script with inputs
srun --cpu-bind=none Rscript $script $num $SLURM_ARRAY_TASK_ID $local
