#downloading all the packages
#if (!require("BiocManager", quietly = TRUE))
# install.packages("BiocManager")

#BiocManager::install("muscle")
library(muscle)

#BiocManager::install("annotate")
library(annotate)

#BiocManager::install("Biostrings")
library(Biostrings)

#installing other packages
#install.packages("seqinr")
library(seqinr)

#install.packages("ape")
library(ape)

#install.packages("dplyr")
library(dplyr)

#install.packages("ggplot2")
library(ggplot2)

#install.packages(reshape2)
library(reshape2)

#install.packages(readr)
library(readr)

#install.packages(stringr)
library(stringr)

#install.packages(littler)
library(littler)

#install.packages(EnvNJ)
library(EnvNJ)

#install.packages(bios2mds)
library(bios2mds)
##########################################################


#muscle reference thing
#creating function to convert fasta to string
fasta_to_string=function(filename){
  fastafile=read.FASTA(filename, type = "DNA")
  DNA_string=fastafile %>% 
    as.character %>% 
    lapply(.,paste0,collapse="") %>%
    unlist %>%
    DNAStringSet
}

#Loading and naming the files
HCov_229E_string=fasta_to_string("./Sequences/HCov-229E_NC_002645")  
HCov_HKU1_string=fasta_to_string("./Sequences/HCov-HKU1_NC_006577") 
HCov_NL63_string=fasta_to_string("./Sequences/HCov-NL63_NC_005831") 
HCov_OC43_string=fasta_to_string("./Sequences/HCov-OC43_NC_006213")
HCov_BS_string=fasta_to_string("./Sequences/hCoV_BS001349.fasta")
HCov_ON_string=fasta_to_string("./Sequences/hCoV_ON078487.fasta")
MERS_string=fasta_to_string("./Sequences/MERS-CoV_NC_019843")
SARS_Cov_string=fasta_to_string("./Sequences/SARS-Cov_NC_004718")
SARS_Cov2_string=fasta_to_string("./Sequences/SARS-CoV2_NC_045512")
BAT_string=fasta_to_string("./Sequences/Bat_Cov_DQ-022305.fasta")
SAMBAR_DEER_string=fasta_to_string("./Sequences/SambarDeerCov_FJ425189.1.fasta")
WATER_DEER_string=fasta_to_string("./Sequences/WaterDeer_MG518518.1.fasta")
WHITE_DEER_OL_string=fasta_to_string("./Sequences/CoV_WTDeer_OL855841.fasta")
WHITE_DEER_FJ_string=fasta_to_string("./Sequences/White-tailedDeerCov_FJ425187.1")
