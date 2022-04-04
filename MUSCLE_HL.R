#Running MUSCLE on sequences

#Installing required packages
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("muscle")
library(BiocManager)
install("genbankr")
install("Biostrings")
install("ggtree")
install("annotate")
install("muscle") 
install("ape")
library(genbankr)
library(ggtree)
library(annotate)
library(ape)
?read.FASTA()
library(muscle)
library(dplyr)
library(Biostrings)
?muscle
#Load Fasta Files as variables
HCoV<- read.FASTA("HCoV-HKU1_NC_006577", type="DNA")
MERS<- read.FASTA("MERS-CoV_NC_019843", type="DNA")
SARS<- read.FASTA("SARS-CoV_NC_004718", type="DNA")

#Convert to strings
HCoVstring<-HCoV %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

MERSstring<-MERS %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

SARSstring<-SARS %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

#Combine sequence strings into single string
DNAstring<-c(HCoVstring, MERSstring, SARSstring)


#Run MUSCLE on sequences
?muscle()
alignment<-muscle::muscle(stringset=DNAstring, quiet=T)
alignment
detail(alignment)

#Plot length on histogram to identify major gaps
seqlen<-as.numeric(lapply(DNAstring, length))
library(ggplot2)
qplot(seqlen)+theme_bw()

#Create distance matrix
alignmentbin <- as.DNAbin(alignment)
AlignDistance<-dist.dna(alignmentbin, model="K80")
class(AlignDistance)
length(AlignDistance)
names(DNAstring)


#View as matrix
AlignDistancemat<-as.matrix(AlignDistance)
dim(AlignDistancemat)
View(AlignDistancemat)
library(reshape2)
DistanceMatrix<-melt(AlignDistancemat)
dim(DistanceMatrix)
View(DistanceMatrix)

#Plot Distance Matrix and change labels
sequences<-c("HCoV-HKU1_NC_006577","MERS-CoV_NC_019843","SARS-CoV_NC_004718")
ggplot(data = DistanceMatrix, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_x_discrete(labels=sequences) + scale_y_discrete(labels=sequences) + scale_fill_gradientn(colours=c("gray87","gray60", "gray50","gray40","gray15"))


