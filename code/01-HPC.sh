#!/bin/bash

# Inputs into bash
num=$1 # Model number to execute

# Fixed inputs
local=0

home=/caldera/hovenweep/projects/usgs/ecosystems/eesc/cscher/iSDM-framework

# Pull species code from MVPv1.csv
while IFS=, read number code model stuff; do
  if [[ $num == $number ]]; then
    spcode=$code
    mod=$model
    break
  fi
done < $home/code/MVPv1.csv

echo "Setting up model: $num $spcode $mod"

# Run model with number input and local=0
sbatch -J $num-$spcode-setup $home/code/01-flexiSDM.sh $num $local


# End script
