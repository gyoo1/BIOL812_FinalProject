# Creating Tree (Neighbour Joining Approach)
CoV_distance<- read.csv(file="CoV_Distance_HL.csv", header=T, row.names=1)
CoV_distance

#Convert to dist object
CoV_distfortree<-as.dist(CoV_distance)

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
