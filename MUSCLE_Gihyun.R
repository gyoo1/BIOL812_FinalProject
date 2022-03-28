library(BiocManager)
library(ape)
library(annotate)
library(muscle)
library(dplyr)

# Read FASTA files
HCov <- read.FASTA("HCoV-HKU1_NC_006577", type = "DNA")
MERS <- read.FASTA("MERS-CoV_NC_019843", type = "DNA")
SARS <- read.FASTA("SARS-CoV_NC_004718", type = "DNA")

HCov
MERS
SARS

# DNASbin to DNAStringSet
HCovDNAstring <- HCov %>% as.character %>% lapply(.,paste0,collapse="") %>% 
        unlist %>% DNAStringSet
MERSDNAstring <- MERS %>% as.character %>% lapply(.,paste0,collapse="") %>% 
        unlist %>% DNAStringSet
SARSDNAstring <- SARS %>% as.character %>% lapply(.,paste0,collapse="") %>% 
        unlist %>% DNAStringSet

HCovDNAstring
MERSDNAstring
SARSDNAstring
