#!/bin/sh

# run findSel.R on bash 

chr=""
group=""
index=0 # change to group index from 0-6
sel="" # type of selection: divergent or balancing
prob=0.99

Rscript findSel.R $chr $group $index $sel $prob
