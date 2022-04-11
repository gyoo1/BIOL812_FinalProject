# Generating Distance Matrix
library(ape)
library(muscle)
library(Biostrings)
alignedreads<- read.FASTA("cov_alignment.fasta", type = "DNA")
AlignDistance<-dist.dna(alignedreads, model="K80")
class(AlignDistance)
length(AlignDistance)

#create distance matrix
distmatrix<- as.matrix(AlignDistance)
distmatrix

#Creating Tree (Neighbour Joining Approach)
#CoV_distance<- read.csv(file="CoV_Distance_HL.csv", header=T, row.names=1) #reading distance matrix
#CoV_distance#might not be necessary with how we have changed our code

#Convert to dist object
CoV_distfortree<-as.dist(distmatrix)

#Create Tree
library(ape)
library(ggtree)
CoVTree<-nj(CoV_distfortree)
str(CoVTree)
class(CoVTree)
ggtree(CoVTree) + geom_tiplab() #
ggtree(CoVTree, layout="radial") + geom_tiplab()

#Save Tree
write.tree(CoVTree, "CoV_NJ_Tree")
