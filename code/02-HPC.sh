#!/bin/bash

# Inputs into bash
num=$1 # Model number to execute

block=$2 # Block to run

days=$3

# Pull species code from MVPv1.csv
while IFS=, read number code model nim stuff; do
  if [[ $num == $number ]]; then
    spcode=$code
    mod=$model
    memnim=$nim
    break
  fi
done < $home/code/MVPv1.csv

echo "Running model: $num $spcode $mod"


# Run model with number input and local=0
sbatch --qos=seven_days_max --mem-per-cpu $memnim --time=$days-00:00:00 --array=1-$chains:1 -J $num-$spcode-$block-nim $home/code/02-flexiSDM.sh $num $block


# End script
