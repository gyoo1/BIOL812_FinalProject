## Distance Matrix and Tree Building ##

# Packages
library(ape)
library(annotate)
library(muscle)
library(dplyr)
library(ggplot2)
library(ggtree)

# Load alignment file
seq_align <- readDNAMultipleAlignment("./Output/cov_alignment.fasta", format = "fasta")

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
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("./Output/DistMatrix.pdf", width = 8, height = 6, units = "in")

# Save Distance Matrix as a .csv file
write.csv(covDMmat, "./Output/DistMatrix.csv")


# Tree Building ----
CovTree <- nj(covDM) #neighbour-joining tree
ggtree(CovTree) #visualize tree
CovTree$tip.label <- sequences #rename tips

#Rectangular Phylogram
ggtree(CovTree, layout="rectangular") + geom_tiplab() + xlim(0, 0.85)

#Circular Cladogram
ggtree(CovTree, layout="circular", branch.length = "none") + 
        xlim(-15, 15) + geom_tiplab(offset = 1)
ggsave("./Output/CoV_cladogram.pdf", width = 30, height = 30, units = "cm") #save tree


# Bootstrapping ----
set.seed(12345)
bs_CovTree <- boot.phylo(nj(covDM), seq_align, FUN = function(x) nj(dist.dna(x)), B = 1000)

#plot bootstrap values on phylogram
ggtree(CovTree, layout="rectangular") + geom_tiplab() + xlim(0, 0.85) +
        geom_text2(aes(subset = !isTip, label = c(1:14, bs_CovTree)), 
                   hjust = 1.1, vjust = -0.6)
ggsave("./Output/CoV_phylogram.pdf", width = 40, height = 20, units = "cm") #save tree

#plot bootstrap values on circular cladogram
ggtree(CovTree, layout="circular", branch.length = "none") + geom_tiplab(offset = 0.2) + 
        geom_point2(aes(subset = (node %in% c(16:26))), shape=16, size=7, 
                    colour = "white") +
        geom_text2(aes(subset = !isTip, label = c(1:14, bs_CovTree)), size = 3)
ggsave("./Output/CoV_cladogram_bs.pdf", width = 30, height = 30, units = "cm") #save tree


# Annotate Tree ----
ggtree(CovTree, branch.length = "none") + geom_text(aes(label=node)) #view node labels

#Get silhouettes from phylopic
species <- c("Odocoileus_virginianus", "Eptesicus_fuscus", "Orthocoronavirinae")
pics <- ggimage::phylopic_uid(species)

#label pics with associated nodes
pics2 <- data.frame(node = c(12, 11, 10, 1, 5, 2), 
                    image = c(rep(pics$uid[1], 4), pics$uid[2], pics$uid[3]),
                    name = c(rep(pics$name[1], 4), pics$name[2], pics$name[3]))

#plot
ggtree(CovTree, layout="rectangular") + geom_tiplab() + xlim(0, 0.85) +
        
        #add boostrap values
        geom_text2(aes(subset = !isTip, label = c(1:14, bs_CovTree)), 
                   hjust = 1.1, vjust = -0.6) +

        #add clade bars
        geom_strip("WtD-CoV", "WD-BCoV", barsize = 2, color = "firebrick2", 
                   label = "Bovine Coronavirus", offset = -0.006, offset.text = 0.01) +
        geom_strip("SARS-CoV", "Bat-SARS-CoV", barsize = 2, color = "springgreen4",
                   label = "SARS", offset = -0.44, offset.text = 0.01) +
        geom_strip("WtD-SARS-CoV2", "SARS-CoV2-BS", barsize = 2, color = "blue2",
                   label = "COVID-19", offset = -0.64, offset.text = 0.01) +
        
        #add pictures
        geom_cladelab(data = pics2[1:3,],
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("grey44"), 
                      offset = 0.05, imagesize = 0.035) +
        geom_cladelab(data = pics2[4,],
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("grey44"), 
                      offset = 0.105, imagesize = 0.04) +
        geom_cladelab(data = pics2[5,],
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("grey44"), 
                      offset = 0.085, imagesize = 0.045) +
        geom_cladelab(data = pics2[6,],
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("blue2"), offset = 0.18)

#save tree
ggsave("./Output/CoV_phylo_annotated.pdf", width = 40, height = 30, units = "cm")


#annotated plot - circular
ggtree(CovTree, layout="circular", branch.length = "none") + geom_tiplab(offset = 0.2) + 
        geom_point2(aes(subset = (node %in% c(16:26))), shape=16, size=7, colour = "white") +
        
        #add bootstrap values
        geom_text2(aes(subset = !isTip, label = c(1:14, bs_CovTree)), size = 3) +
        
        #add clade bars
        geom_strip("WtD-CoV", "WD-BCoV", barsize = 2, color = "firebrick2", 
                   label = "Bovine Coronavirus", offset = 7.9, offset.text = 0.8) +
        geom_strip("SARS-CoV", "Bat-SARS-CoV", barsize = 2, color = "springgreen4",
                   label = "SARS", offset = 9.1, offset.text = 3) +
        geom_strip("WtD-SARS-CoV2", "SARS-CoV2-BS", barsize = 2, color = "blue2",
                   label = "COVID-19", offset = 11, offset.text = 2) +
        
        #add pictures
        geom_cladelab(data = pics2[1:3,],
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("grey44"), 
                      offset = 5.45, imagesize = 0.035) +
        geom_cladelab(data = pics2[4,],
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("grey44"), 
                      offset = 8.8, imagesize = 0.04) +
        geom_cladelab(data = pics2[5,],
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("grey44"), 
                      offset = 7, imagesize = 0.045) +
        geom_cladelab(data = pics2[6,],
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("blue2"), offset = 13)

#save tree
ggsave("./Output/CoV_clado_annotated.pdf", width = 30, height = 30, units = "cm")
