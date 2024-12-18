#!/bin/bash
cd cmake-build-debug/
program=./Pythia-project-tree-creator
num_terminals=$(grep -c ^processor /proc/cpuinfo)
for ((i=1; i<=num_terminals; i++)); do
    x_offset=$((i * 50))
    y_offset=0
    xterm -geometry 50x50+${x_offset}+${y_offset} -bg black -fg green -e "$program" &
done
wait
cd ../results/
rm d0.root 
hadd d0.root *.root
