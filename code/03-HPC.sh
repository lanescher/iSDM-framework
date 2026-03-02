#!/bin/bash

# Inputs into bash
num=$1 # Model number to execute

block=$2

# Pull species code from MVPv1.csv
while IFS=, read number code model nim sum stuff; do
  if [[ $num == $number ]]; then
    spcode=$code
    mod=$model
    memsum=$sum
    break
  fi
done < code/MVPv1.csv

echo "Summarizing model: $num $spcode $mod"
echo "Block: $block"



sbatch --qos=seven_days_max --mem-per-cpu $memsum --time=0-12:00:00 -J $num-$spcode-$block-sum code/03-flexiSDM.sh $num $block

# End script
