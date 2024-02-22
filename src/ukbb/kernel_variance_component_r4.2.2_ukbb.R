#!/usr/bin/env Rscript

## read arguments
args=(commandArgs(TRUE))
gdsfile=as.character(args[1])
groupfile=as.character(args[2])
phenfile=as.character(args[3])
ID_col=as.character(args[4])
nullfile=as.character(args[5])
outfile=as.character(args[6])

## specify the packages
.libPaths(c(“rpackages4_2_2”,.libPaths()))

## load source scripts
source(“UKBB_200KWES_CVD/GENESIS_adaptation_source.R”)
source(“meta_transcript/src/ExtractKernalStatistics_SPA_transcript_UKBB.R”)

## perfrom analysis
kernell_variance_component_ukbb(gdsfile=gdsfile,groupfile=groupfile,phenfile=phenfile,ID_col=ID_col,nullfile=nullfile,outfile=outfile, test=“ExtractKernelStatistics”, vc.test=“Score.SPA”, AF.max=0.001, MAC.max=Inf,use.weights=FALSE)

##quit!
sessionInfo()
quit("no")