#!/usr/bin/env Rscript

## read arguments
args=(commandArgs(TRUE))
gdsfile=as.character(args[1])
varfile=as.character(args[2])
phenfile=as.character(args[3])
ID_col=as.character(args[4])
nullfile=as.character(args[5])
test=as.character(args[6])
outfile=as.character(args[7])


## specify the packages
.libPaths(c("rpackages4_2_2",.libPaths()))

## load source scripts
source("UKBB_200KWES_CVD/GENESIS_adaptation_source.R")
source("meta_transcript/src/singlevariant_r4.2.2_ukbb.R")

## perfrom analysis
singleassoc_varinfo_ukbb(gdsfile=gdsfile,varfile=varfile,phenfile=phenfile,ID_col="scanID",nullfile=nullfile,stest=test,outfile=outfile)

##quit!
sessionInfo()
quit("no")