install.packages("BiocManager")
install.packages("ape")
install.packages("dplyr")

library(BiocManager)
library(ape)
library(dplyr)

#read sequences 
HCov <- read.FASTA("./Sequences/HCoV-HKU1_NC_006577", type="DNA")
MERS <- read.FASTA("./Sequences/MERS-CoV_NC_019843", type="DNA")
SARS <- read.FASTA("./Sequences/SARS-CoV_NC_004718", type="DNA")
HCovNL63 <- read.FASTA("./Sequences/HCov-NL63_NC_005831", type="DNA")
HCovOC43 <- read.FASTA("./Sequences/HCov-OC43_NC_006213", type="DNA")
HCov2 <- read.FASTA("./Sequences/SARS-CoV2_NC_045512", type="DNA")
HCov229E <- read.FASTA("./Sequences/HCoV-229E_NC_002645", type="DNA")
Bat <- read.FASTA("./Sequences/Bat_Cov_DQ-022305.fasta", type="DNA")
Deer <- read.FASTA("./Sequences/CoV_WTDeer_OL855841.fasta", type="DNA")
Sambar <- read.FASTA("./Sequences/SambarDeerCov_FJ425189.1.fasta", type="DNA")
Wdeer <- read.FASTA("./Sequences/WaterDeer_MG518518.1.fasta", type="DNA")
Whdeer <- read.FASTA("./Sequences/White-tailedDeerCov_FJ425187.1", type="DNA")
HCovBS <- read.FASTA("./Sequences/hCoV_BS001349.fasta", type="DNA")
HCovON <- read.FASTA("./Sequences/hCoV_ON078487.fasta", type="DNA")

#create string set
HCovstring <- HCov %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
MERSstring <- MERS %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
SARSstring <- SARS %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
HCovNL63string <- HCovNL63 %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
HCovOC43string <- HCovOC43 %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
HCov2string <- HCov2 %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
HCov229Estring <- HCov229E %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
Batstring <- Bat %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
Deerstring <- Deer %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
Sambarstring <- Sambar %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
Wdeerstring <- Wdeer %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
Whdeerstring <- Whdeer %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
HCovBSstring <- HCovBS %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
HCovONstring <- HCovON %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet

HCovstring
MERSstring 
SARSstring
HCovNL63string 
HCovOC43string 
HCov2string 
HCov229Estring
Batstring
Deerstring
Sambarstring
Wdeerstring 
Whdeerstring 
HCovBSstring 
HCovONstring 

#muscle 
library(muscle)
DNAset <- c(HCovstring,MERSstring,SARSstring,HCovNL63string,HCovOC43string,HCov2string,HCov229Estring,Batstring,Deerstring,Sambarstring,Wdeerstring,Whdeerstring,HCovBSstring,HCovONstring)
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



