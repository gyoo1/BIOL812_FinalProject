install.packages("BiocManager")
install.packages("ape")
install.packages("dplyr")

library(BiocManager)
library(ape)
library(dplyr)

#read sequences 
HCov <- read.FASTA("HCoV-HKU1_NC_006577", type="DNA")
MERS <- read.FASTA("MERS-CoV_NC_019843", type="DNA")
SARS <- read.FASTA("SARS-CoV_NC_004718", type="DNA")

#create string set
HCovstring <- HCov %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
MERSstring <- MERS %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
SARSstring <- SARS %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet

HCovstring
MERSstring 
SARSstring

#muscle 
library(muscle)
DNAset <- c(HCovstring,MERSstring,SARSstring)
DNAalign <- muscle(stringset = DNAset, quiet = T)
DNAalign

#see alignment 
detail(DNAalign)

#Histogram of DNAset
Seq <-as.numeric(lapply(DNAset,length))
library(ggplot2)
qplot(Seq)+theme_bw()

#distance matrix
align <- as.DNAbin(DNAalign)
DNAsetDM<-dist.dna(align, model="K80")
class(DNAsetDM)

DNADM<-as.matrix(DNAsetDM)
dim(DNADM)

library(reshape2)
DM<-melt(DNADM)
dim(DM)

#plot
ggplot(data = DM, aes(x=Var1, y=Var2, fill=value)) + 
  geom_tile() + scale_fill_gradientn(colours=c("black","grey1","white")) +
  scale_x_discrete(labels = c("HCoV-HKU1_NC_006577","MERS-CoV_NC_019843","SARS-CoV_NC_004718")) + scale_y_discrete(labels = c("HCoV-HKU1_NC_006577","MERS-CoV_NC_019843","SARS-CoV_NC_004718")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))


