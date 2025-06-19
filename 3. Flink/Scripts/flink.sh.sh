#!/bin/bash

dir_flink="/home/aliciaportela/AntropoGeo/flink"

$dir_flink/flink task=estimate numThreads=10 B="(-2.0,1.8)" A_max="4.0" lnMu_g="(-4.0,-0.0)" lnNu_g="(-5.0,-0.0)" s_max=14 beta="(-2.0,1.8)" alpha_max=4.0 numIterations=500000 burnin=300000 thinning=100 lnkappa="(-10.0,-0.1)" logFile=logfile_human sigmaProp_mu=0.005 sigmaProp_nu=0.05 sigmaProp_kappa=0.05 data=pop_chr12_p.txt

