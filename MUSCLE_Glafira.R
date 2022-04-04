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

HCov_229E=read.FASTA("./Sequences/HCov-229E_NC_002645", type = "DNA")
HCov_HKU1=read.FASTA("./Sequences/HCov-HKU1_NC_006577", type = "DNA")
HCov_NL63=read.FASTA("./Sequences/HCov-NL63_NC_005831", type = "DNA")
HCov_OC43=read.FASTA("./Sequences/HCov-OC43_NC_006213", type = "DNA")
HCov_BS=read.FASTA("./Sequences/hCoV_BS001349.fasta", type = "DNA")
HCov_ON=read.FASTA("./Sequences/hCoV_ON078487.fasta", type = "DNA")
MERS=read.FASTA("./Sequences/MERS-CoV_NC_019843", type = "DNA")
SARS_Cov=read.FASTA("./Sequences/SARS-Cov_NC_004718", type = "DNA")
SARS_Cov2=read.FASTA("./Sequences/SARS-CoV2_NC_045512", type = "DNA")
BAT=read.FASTA("./Sequences/Bat_Cov_DQ-022305.fasta", type = "DNA")
SAMBAR_DEER=read.FASTA("./Sequences/SambarDeerCov_FJ425189.1.fasta", type = "DNA")
Water_DEER=read.FASTA("./Sequences/WaterDeer_MG518518.1.fasta", type = "DNA")
WHITE_DEER_OL=read.FASTA("./Sequences/CoV_WTDeer_OL855841.fasta", type = "DNA")
WHITE_DEER_FJ=read.FASTA("./Sequences/White-tailedDeerCov_FJ425187.1", type = "DNA")

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

#saving the matrix to make the tree 
write.csv(Alignment_mat,file="cov_distance_GE.csv")
