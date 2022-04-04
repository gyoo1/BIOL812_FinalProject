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
HCovstring <- HCov %>% as.character %>% lapply(.,paste0,collapse="") %>%
        unlist %>% DNAStringSet
MERSstring <- MERS %>% as.character %>% lapply(.,paste0,collapse="") %>%
        unlist %>% DNAStringSet
SARSstring <- SARS %>% as.character %>% lapply(.,paste0,collapse="") %>%
        unlist %>% DNAStringSet

HCovstring
MERSstring
SARSstring

# Multiple alignments with MUSCLE
DNAstring <- c(HCovstring, MERSstring, SARSstring)
COVAlign <- muscle::muscle(stringset= DNAstring, quiet=T)
COVAlign

# Identify Gaps
SeqLen<-as.numeric(lapply(DNAstring,length))
library(ggplot2)
qplot(SeqLen)+theme_bw() #no large gaps

# Distance Matrix ----

# Convert to DNAbin
COVAlign <- as.DNAbin(COVAlign)
COVDM <- dist.dna(COVAlign, model="K80")
class(COVDM)
length(COVDM)

# Convert to linear matrix
library(reshape2)
COVDMmat<-as.matrix(COVDM)
PDat <- melt(COVDMmat)
dim(PDat)
View(PDat)

sequences <- c("HCoV-HKU1_NC_006577", "MERS-CoV_NC_019843", "SARS-CoV_NC_004718")
ggplot(data = PDat, aes(x=Var1, y=Var2, fill=value)) + geom_tile() +
        scale_x_discrete(labels = sequences) + scale_y_discrete(labels = sequences) +
        xlab("Strain") + ylab("Strain") + theme_bw() + 
        scale_fill_gradientn(colours = c("white", "white", "white", "white",
                                         "white", "white", "white", "white",
                                         "white", "white", "white", "green"))
