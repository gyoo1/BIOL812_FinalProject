install.packages("BiocManager")
install.packages("ape")
install.packages("dplyr")


library(BiocManager)
library(ape)
library(muscle)
library(dplyr)

HCov <- read.FASTA("HCoV-HKU1_NC_006577", type="DNA")
MERS <- read.FASTA("MERS-CoV_NC_019843", type="DNA")
SARS <- read.FASTA("SARS-CoV_NC_004718", type="DNA")
