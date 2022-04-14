## NJ Tree for SARS-CoV-2 Strains Only ##

# Packages
library(ape)
library(EnvNJ)
library(Biostrings)
library(dplyr)
library(annotate)
library(muscle)
library(bios2mds)
library(msa)
library(ggplot2)
library(ggtree)
library(ggimage)

# Read FASTA files ----
CoV2 <- c("CoV_WTDeer_OL855841", "hCoV_BS001349", "hCoV_ON078487", "SARS-CoV2_NC_045512")
fastaconc(CoV2, inputdir = "./Sequences", out.file = "./Output/CoV2_concatenated.fasta") 
#concatenate all sequences into a single .fasta file

CoV2_seq <- read.FASTA("./Output/CoV2_concatenated.fasta") #read concatenated sequence

# Convert DNASbin to DNAStringSet file ----
CoV2_string <- CoV2_seq %>% as.character %>% lapply(.,paste0,collapse="") %>%
        unlist %>% DNAStringSet

# Multiple alignments with MUSCLE ----
CoV2_align <- muscle::muscle(stringset= CoV2_string, quiet=T)

# Save alignment output as fasta file ----
CoV2_align_as_align <- msaConvert(CoV2_align, "bios2mds::align") #convert to align object
export.fasta(CoV2_align_as_align, outfile = "./Output/CoV2_alignment.fasta")


# Distance Matrix ----

# Convert to DNAbin
CoV2_align <- as.DNAbin(CoV2_align)

# Distance Matrix
cov2DM <- dist.dna(CoV2_align, model="K80")
class(cov2DM) #check class
length(cov2DM) #check length

# Convert to linear matrix
library(reshape2)
cov2DMmat <- as.matrix(cov2DM)
PDat2 <- melt(cov2DMmat)
dim(PDat2) #check dimensions

# Plot Distance Matrix
CoV2_sequences <- c("WtD-SARS-CoV2", "SARS-CoV2-BS", "SARS-CoV2-Omicron", "SARS-CoV2-NC") 
#rename sequences

#plot
ggplot(data = PDat2, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
        scale_x_discrete(labels = sequences) + scale_y_discrete(labels = sequences) +
        xlab("Strain") + ylab("Strain") + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
ggsave("./Output/DistMatrix_CoV2.pdf", width = 8, height = 6, units = "in")

# Save Distance Matrix as a .csv file
write.csv(cov2DMmat, "DistMatrix_CoV2.csv")


# Tree Building ----
Cov2Tree <- nj(cov2DM) #neighbour-joining tree
ggtree(Cov2Tree) #visualize tree
Cov2Tree$tip.label <- CoV2_sequences #rename tips

#Rectangular Phylogram
ggtree(Cov2Tree, layout="rectangular") + geom_tiplab() + xlim(NA, 0.0015)

#Circular Cladogram
ggtree(Cov2Tree, layout="circular", branch.length = "none") + 
        xlim(-15, 15) + geom_tiplab(offset = 1)

# Bootstrapping ----
set.seed(12345)
bs_Cov2Tree <- boot.phylo(nj(cov2DM), CoV2_align, FUN = function(x) nj(dist.dna(x)), B = 1000)

#plot bootstrap values on cladogram ----
species1 <- "Orthocoronavirinae"
pics3 <- ggimage::phylopic_uid(species1)
pics4 <- data.frame(node = 5, image = pics3$uid, name = pics3$name) 
        #label pics with associated nodes

#plot
ggtree(Cov2Tree, layout="circular", branch.length = "none") + 
        xlim(-15, 15) + geom_tiplab(offset = 1) +
        geom_text2(aes(subset = !isTip, label = c(1:4, bs_Cov2Tree)),
                   nudge_x = -1, nudge_y = 0.08) +
        geom_cladelab(data = pics4,
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("blue2"),
                      imagesize = 0.25,
                      offset = -17)
        
ggsave("./Output/CoV2_cladogram_annotated.pdf", width = 20, height = 20, units = "cm") #save tree


# Facets ----
#run Tree_Gihyun.R first
#requires plot from Tree_Gihyun.R loaded into the environment
FullTree <- ggtree(CovTree, layout="circular", branch.length = "none") + geom_tiplab(offset = 0.2) + 
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

SubTree <- ggtree(Cov2Tree, layout="circular", branch.length = "none") + 
        xlim(-15, 15) + geom_tiplab(offset = 1) +
        geom_text2(aes(subset = !isTip, label = c(1:4, bs_Cov2Tree)),
                   nudge_x = -1, nudge_y = 0.08) +
        geom_cladelab(data = pics4,
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("blue2"),
                      imagesize = 0.25,
                      offset = -17)

library(grid)
pdf("./Output/Covid_Cladograms.pdf", width = 18, height = 12)
        grid.newpage() # Open new page on grid device
        pushViewport(viewport(layout = grid.layout(3, 5))) #set up plotting grid with 1 row + 2 cols
        print(FullTree, vp = viewport(layout.pos.row = 1:3,layout.pos.col = 1:3))
        print(SubTree, vp = viewport(layout.pos.row = 1:3,layout.pos.col = 4:5))
dev.off()

