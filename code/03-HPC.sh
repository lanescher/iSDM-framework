#!/bin/bash

# Inputs into bash
num=$1 # Model number to execute

block=$2

# Fixed inputs
local=0
maxchain=3


# Set working directory
home=/caldera/hovenweep/projects/usgs/ecosystems/eesc/cscher/iSDM-framework


# Pull species code from MVPv1.csv
while IFS=, read number code model nim sum stuff; do
  if [[ $num == $number ]]; then
    spcode=$code
    mod=$model
    memsum=$sum
    break
  fi
done < $home/code/MVPv1.csv

echo "Summarizing model: $num $spcode $mod"
echo "Block: $block"
echo "Local: $local"



sbatch --qos=seven_days_max --mem-per-cpu $memsum --time=0-12:00:00 -J $num-$spcode-$block-sum $home/code/03-flexiSDM.sh $num $block $maxchain $local

# End script
