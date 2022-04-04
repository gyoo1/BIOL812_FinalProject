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
HCoV<-read.FASTA("./Sequences/HCoV-HKU1_NC_006577", type="DNA")
MERS<-read.FASTA("./Sequences/MERS-CoV_NC_019843", type="DNA")
SARS<-read.FASTA("./Sequences/SARS-CoV_NC_004718", type="DNA")
BAT<- read.FASTA("./Sequences/Bat_Cov_DQ-022305.fasta", type="DNA")
DEER_OL<- read.FASTA("./Sequences/CoV_WTDeer_OL855841.fasta", type="DNA")
HCoV_BS<- read.FASTA("./Sequences/hCoV_BS001349.fasta", type="DNA")
HCoV_ON<- read.FASTA("./Sequences/hCoV_ON078487.fasta", type="DNA")
HCoV_229E<- read.FASTA("./Sequences/HCoV-229E_NC_002645", type="DNA")
HCoV_NL <- read.FASTA("./Sequences/HCov-NL63_NC_005831", type="DNA")
HCoV_OC<- read.FASTA("./Sequences/HCov-OC43_NC_006213", type="DNA")
SAMBAR<- read.FASTA("./Sequences/SambarDeerCov_FJ425189.1.fasta", type="DNA")
SARS_CoV2<- read.FASTA("./Sequences/SARS-CoV2_NC_045512", type="DNA")
DEER_MG<- read.FASTA("./Sequences/WaterDeer_MG518518.1.fasta", type="DNA")
DEER_FJ<- read.FASTA("./Sequences/White-tailedDeerCov_FJ425187.1", type="DNA")

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

BATstring<-BAT %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

DEER_OLstring<-DEER_OL %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

HCoV_BSstring<-HCoV_BS %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

HCoV_ONstring<-HCoV_ON %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

HCoV_229Estring<-HCoV_229E %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

HCoV_NLstring<-HCoV_NL %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

HCoV_OCstring<-HCoV_OC %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

SAMBARstring<-SAMBAR %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

SARS_C0V2string<-SARS_CoV2 %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

DEER_MGstring<-DEER_MG %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

DEER_FJstring<-DEER_FJ %>%
  as.character %>%
  lapply(.,paste0,collapse="") %>%
  unlist%>%
  DNAStringSet

#Combine sequence strings into single string
DNAstring<-c(HCoVstring, MERSstring, SARSstring, BATstring, DEER_OLstring,
             HCoV_BSstring, HCoV_ONstring, HCoV_229Estring, HCoV_NLstring, 
             HCoV_OCstring, SAMBARstring, SARS_C0V2string, DEER_MGstring, DEER_FJstring)


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

#Save Distance Matrix Values in Local File
write.csv(AlignDistancemat, file="CoV_Distance_HL.csv")


