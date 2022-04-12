## Distance Matrix and Tree Building ##

# Packages
library(BiocManager)
library(ape)
library(annotate)
library(muscle)
library(dplyr)
library(ggplot2)
library(ggtree)

# Load alignment file
seq_align <- readDNAMultipleAlignment("cov_alignment.fasta", format = "fasta")

# Distance Matrix ----

# Convert to DNAbin
seq_align <- as.DNAbin(seq_align)

# Distance Matrix
covDM <- dist.dna(seq_align, model="K80")
class(covDM) #check class
length(covDM) #check length

# Convert to linear matrix
library(reshape2)
covDMmat <- as.matrix(covDM)
PDat <- melt(covDMmat)
dim(PDat) #check dimensions

# Plot Distance Matrix
sequences <- c("WtD-SARS-CoV2", "SARS-CoV2-Omicron", "SARS-CoV2-BS", "SARS-CoV2-NC", 
               "Bat-SARS-CoV", "SARS-CoV", "MERS-CoV", "HCoV-HKU1", "HCoV-OC43", "WD-BCoV",
               "SD-CoV", "WtD-CoV", "HCoV-NL63", "HCoV-229E") #rename sequences
#plot
ggplot(data = PDat, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
        scale_x_discrete(labels = sequences) + scale_y_discrete(labels = sequences) +
        xlab("Strain") + ylab("Strain") + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_fill_gradientn(colours = c("azure", "cyan3", "royalblue3"))

# Save Distance Matrix as a .csv file
write.csv(covDMmat, "Cov_Distance.csv")


# Tree Building ----
CovTree <- nj(covDM) #neighbour-joining tree
ggtree(CovTree) #visualize tree
CovTree$tip.label <- sequences #rename tips

#Rectangular Phylogram
ggtree(CovTree, layout="rectangular") + geom_tiplab() + xlim(0, 0.85)

#Circular Cladogram
ggtree(CovTree, layout="circular", branch.length = "none") + 
        xlim(-15, 15) + geom_tiplab(offset = 1)
ggsave("CoV_cladogram.pdf", width = 30, height = 30, units = "cm") #save cladogram


# Bootstrapping ----
bs_CovTree <- boot.phylo(nj(covDM), seq_align, FUN = function(x) nj(dist.dna(x)))

#plot bootstrap values on phylogram
ggtree(CovTree, layout="rectangular") + geom_tiplab() + xlim(0, 0.85) +
        geom_text2(aes(subset = !isTip, label = c(1:14, bs_CovTree)), hjust = 1.1, vjust = 1.2)
ggsave("CoV_phylogeny.pdf", width = 40, height = 20, units = "cm") #save phylogram

#plot bootstrap values on cladogram
ggtree(CovTree, layout="circular", branch.length = "none") + geom_tiplab(offset = 0.2) + 
        geom_text2(aes(subset = !isTip, label = c(1:14, bs_CovTree)),
                   nudge_x = -0.45, nudge_y = -0.19)
ggsave("CoV_cladogram_boot.pdf", width = 30, height = 30, units = "cm") #save cladogram


# Annotate Tree ----
ggtree(CovTree, branch.length = "none") + geom_text(aes(label=node)) #view node labels

ggtree(CovTree) + xlim(0, 0.85) + geom_tiplab() +
        geom_strip("WtD-CoV", "WD-BCoV", barsize = 2, color = "red", 
                   label = "Bovine Coronavirus", offset = -0.05, offset.text = 0.01)
ggsave("CoV_phylogeny_annotated.pdf", width = 40, height = 20, units = "cm") #save cladogram

