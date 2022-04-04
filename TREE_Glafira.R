#Using neighbour joining to make a tree

#downloading all the packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("muscle")
library(muscle)

BiocManager::install("annotate")
library(annotate)

BiocManager::install("Biostrings")
library(Biostrings)

#installing other packages
install.packages("seqinr")
library(seqinr)

install.packages("ape")
library(ape)

install.packages("dplyr")
library(dplyr)

install.packages("ggplot2")
library(ggplot2)

install.packages("reshape2")
library(reshape2)


#uploading the distance matrix values made in the MUSCLE script 
Alignment_mat_fortree=read.csv("cov_distance_GE.csv",header=T, row.names=1)
View(Alignment_mat_fortree)
CovDist_fortree=as.dist(Alignment_mat_fortree)

#ok! lets make the tree
library(BiocManager)
install("ggtree")
library(ggtree)
Cov_Tree<-nj(CovDist_fortree)
str(Cov_Tree)
class(Cov_Tree)
ggtree(Cov_Tree,layout="circular") + geom_tiplab()
