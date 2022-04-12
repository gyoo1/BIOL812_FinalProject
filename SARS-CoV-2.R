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
filenames <- c("CoV_WTDeer_OL855841", "hCoV_BS001349", "hCoV_ON078487", "SARS-CoV2_NC_045512")
fastaconc(filenames, inputdir = "./Sequences", out.file = "./CoV2_concatenated.fasta") 
#concatenate all sequences into a single .fasta file

seq <- read.FASTA("CoV2_concatenated.fasta") #read concatenated sequence

# Convert DNASbin to DNAStringSet file ----
seq_string <- seq %>% as.character %>% lapply(.,paste0,collapse="") %>%
        unlist %>% DNAStringSet

# Multiple alignments with MUSCLE ----
seq_align <- muscle::muscle(stringset= seq_string, quiet=T)

# Save alignment output as fasta file ----
seq_align_as_align <- msaConvert(seq_align, "bios2mds::align") #convert to align object
export.fasta(seq_align_as_align, outfile = "CoV2_alignment.fasta")


# Distance Matrix ----

# Convert to DNAbin
seq_align <- as.DNAbin(seq_align)

# Distance Matrix
cov2DM <- dist.dna(seq_align, model="K80")
class(cov2DM) #check class
length(cov2DM) #check length

# Convert to linear matrix
library(reshape2)
cov2DMmat <- as.matrix(cov2DM)
PDat <- melt(cov2DMmat)
dim(PDat2) #check dimensions

# Plot Distance Matrix
sequences <- c("WtD-SARS-CoV2", "SARS-CoV2-BS", "SARS-CoV2-Omicron", "SARS-CoV2-NC") 
#rename sequences

#plot
ggplot(data = PDat2, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
        scale_x_discrete(labels = sequences) + scale_y_discrete(labels = sequences) +
        xlab("Strain") + ylab("Strain") + theme_bw() + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
        scale_fill_gradientn(colours = c("azure", "cyan3", "royalblue3"))

# Save Distance Matrix as a .csv file
write.csv(cov2DMmat, "Cov2_Distance.csv")


# Tree Building ----
Cov2Tree <- nj(cov2DM) #neighbour-joining tree
ggtree(Cov2Tree) #visualize tree
Cov2Tree$tip.label <- sequences #rename tips

#Rectangular Phylogram
ggtree(Cov2Tree, layout="rectangular") + geom_tiplab() + xlim(NA, 0.0015)

#Circular Cladogram
ggtree(Cov2Tree, layout="circular", branch.length = "none") + 
        xlim(-15, 15) + geom_tiplab(offset = 1)

# Bootstrapping ----
bs_Cov2Tree <- boot.phylo(nj(cov2DM), seq_align, FUN = function(x) nj(dist.dna(x)))

#plot bootstrap values on phylogram
ggtree(Cov2Tree, layout="rectangular") + geom_tiplab() + xlim(NA, 0.0015) +
        geom_text2(aes(subset = !isTip, label = c(1:4, bs_Cov2Tree)), hjust = 1.1, vjust = 1.2)
ggsave("CoV2_phylogeny.pdf", width = 40, height = 20, units = "cm") #save phylogram

#plot bootstrap values on cladogram
species <- "Orthocoronavirinae"
pics <- ggimage::phylopic_uid(species)
pics2 <- data.frame(node = 5, image = pics$uid, name = pics$name) 
        #label pics with associated nodes

#plot
ggtree(Cov2Tree, layout="circular", branch.length = "none") + 
        xlim(-15, 15) + geom_tiplab(offset = 1) +
        geom_text2(aes(subset = !isTip, label = c(1:4, bs_Cov2Tree)),
                   nudge_x = -1, nudge_y = 0.08) +
        geom_cladelab(data = pics2,
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("blue"),
                      imagesize = 0.25,
                      offset = -17)
        
ggsave("CoV2_cladogram.pdf", width = 30, height = 30, units = "cm") #save cladogram


# Facets ----
#run Tree_Gihyun.R first
#requires plot from Tree_Gihyun.R loaded into the environment
FullTree <- ggtree(CovTree, layout="circular", branch.length = "none") + #all strains
        geom_tiplab(offset = 0.2) + 
        geom_text2(aes(subset = !isTip, label = c(1:14, bs_CovTree)),
                   nudge_x = -0.45, nudge_y = -0.19)
SubTree <- ggtree(Cov2Tree, layout="circular", branch.length = "none") + #only SARS-CoV-2
        xlim(-15, 15) + geom_tiplab(offset = 1) +
        geom_text2(aes(subset = !isTip, label = c(1:4, bs_Cov2Tree)),
                   nudge_x = -1, nudge_y = 0.08) +
        geom_cladelab(data = pics2,
                      mapping = aes(node = node, label = name, image = image), 
                      geom = "phylopic", imagecolor = c("blue"),
                      imagesize = 0.25,
                      offset = -17)

library(grid)
grid.newpage() # Open new page on grid device
pushViewport(viewport(layout = grid.layout(1, 2))) #set up plotting grid with 1 row + 2 cols
print(FullTree, vp = viewport(layout.pos.row = 1,layout.pos.col = 1))
print(SubTree, vp = viewport(layout.pos.row = 1,layout.pos.col = 2))

