#!/usr/bin/env Rscript

library(seqinr)

source("/home/tpillone/Dropbox/phylogenomics/bin/getPttV1.2.R")
source("/home/tpillone/Dropbox/phylogenomics/bin/getGenesV1.10.R")
source("/home/tpillone/Dropbox/phylogenomics/bin/detectionV3.2.R")

args = commandArgs(TRUE)

ptt<-getPtt(args[1])
fna<-read.fasta(args[2])
genes<-getGenes(fna,ptt,stat=TRUE)

write.table(genes,file="",sep="\t",col.names = F,row.names=F)
