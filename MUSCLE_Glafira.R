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

install.packages(reshape2)
library(reshape2)

#Loading and naming the files
HCov=read.FASTA("HCov-HKU1_NC_006577", type = "DNA")
MERS=read.FASTA("MERS-CoV_NC_019843", type = "DNA")
SARS=read.FASTA("SARS-Cov_NC_004718", type = "DNA")

#converting to strings

#HCovstring
HCovstring=HCov %>% 
  as.character %>% 
  lapply(.,paste0,collapse="") %>%
  unlist %>%
  DNAStringSet
head(HCovstring)
  
#MERSstring

MERSstring=MERS %>% 
  as.character %>% 
  lapply(.,paste0,collapse="") %>%
  unlist %>%
  DNAStringSet
head(MERSstring)

#MERSstring

SARSstring=SARS %>% 
  as.character %>% 
  lapply(.,paste0,collapse="") %>%
  unlist %>%
  DNAStringSet
head(SARSstring)

#Using muscle to align the sequences
stringset=c(HCovstring,MERSstring,SARSstring)
Alignment<-muscle::muscle(stringset=c(HCovstring,MERSstring,SARSstring), quiet=T)
Alignment

#check the gaps in alignment sequence 
Seqlength<-as.numeric(lapply(stringset,length))
qplot(Seqlength)+theme_bw()

#looks like the alignment doesn't have a lot of gaps, so we will not remove them
#doing the distance matrix
AlignmentBin <- as.DNAbin(Alignment)
AlignmentBin_DM<-dist.dna(AlignmentBin, model="K80")
class(AlignmentBin_DM)
length(AlignmentBin_DM)
Alignment_mat<-as.matrix(AlignmentBin_DM)
dim(Alignment_mat)
View(Alignment_mat)
library(reshape2)
Alignment_P<-melt(Alignment_mat)
dim(Alignment_P)
View(Alignment_P)

#plotting the matrix
ggplot(data = Alignment_P, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile()+
  scale_x_discrete(labels=c("HCov","MERS","SARS"))+
  scale_y_discrete(labels=c("HCov","MERS","SARS"))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))+
  scale_fill_gradientn(colours=c("white","grey96","grey88","grey76","grey77","grey17"))
  