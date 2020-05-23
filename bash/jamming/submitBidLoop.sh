#!/bin/bash

# loop parameters
calA0start=1.04
calA0step=0.02
calA0end=1.18

kb=0

partition=scavenge

for calA0 in `seq $calA0start $calA0step $calA0end`; do
    bash bidCellJamming.sh 16 16 $calA0 1.0 $kb $partition 0-18:00:00 1 100 1
done