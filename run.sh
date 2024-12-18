#!/bin/bash
cd cmake-build-debug/
program=./Pythia-project-tree-creator
num_terminals=$(grep -c ^processor /proc/cpuinfo)
for ((i=1; i<=num_terminals; i++)); do
    x_offset=0
    y_offset=$((i * 20))
    xterm -geometry 50x20+${x_offset}+${y_offset} -bg black -fg green -hold -e "$program" &
done

