## Distance Matrix and Tree Building ##

# Packages
library(BiocManager)
library(ape)
library(annotate)
library(muscle)
library(dplyr)
library(ggtree)

# Load alignment file
seq_align <- readDNAMultipleAlignment("cov_alignment.fasta", format = "fasta")

# Distance Matrix ----

# Convert to DNAbin
seq_align <- as.DNAbin(seq_align)
covDM <- dist.dna(seq_align, model="K80")
class(covDM)
length(covDM)

# Convert to linear matrix
library(reshape2)
covDMmat <- as.matrix(covDM)
PDat <- melt(covDMmat)
dim(PDat)
View(PDat)

sequences <- c("HCoV-HKU1_NC_006577", "MERS-CoV_NC_019843", "SARS-CoV_NC_004718")
ggplot(data = PDat, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
        scale_x_discrete(labels = sequences) + scale_y_discrete(labels = sequences) +
        xlab("Strain") + ylab("Strain") + theme_bw() + 
        scale_fill_gradientn(colours = c("white", "white", "white", "white",
                                         "white", "white", "white", "white",
                                         "white", "white", "white", "green"))


# Save Distance Matrix as .csv
write.csv(covDMmat, "Cov_Distance.csv")

# Read in Distance
CovDistMat <- read.csv("Cov_Distance.csv", header = T, row.names = 1)
CovDist <- as.dist(CovDistMat)

# Tree
CovTree <- nj(CovDist)
ggtree(CovTree)

ggtree(CovTree,layout="rectangular") + geom_tiplab() #label tips
ggtree(CovTree,layout="circular") + geom_tiplab() #circular phylogeny



