#!/bin/sh

# run findSel.R on bash 
# run it for each group (and hierarchy) and each type of selection

chr="chr12" 
group="" # name of the group
index=0 # change to group index from 0-6
sel="" # type of selection: divergent or balancing
prob=0.99999

Rscript findSel.R $chr $group $index $sel $prob
