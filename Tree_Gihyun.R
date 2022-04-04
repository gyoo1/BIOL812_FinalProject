library(BiocManager)
library(ape)
library(annotate)
library(muscle)
library(dplyr)
library(ggtree)

# Read in Distance
CovDistMat <- read.csv("Cov_Distance.csv", header = T, row.names = 1)
CovDist <- as.dist(CovDistMat)

# Tree
CovTree <- nj(CovDist)
ggtree(CovTree)

ggtree(CovTree,layout="rectangular") + geom_tiplab() #label tips
ggtree(CovTree,layout="circular") + geom_tiplab() #circular phylogeny



