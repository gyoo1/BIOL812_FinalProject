library(BiocManager)
library(ape)
library(annotate)
library(muscle)
library(dplyr)

# Read FASTA files
Bat <- read.FASTA("./Sequences/Bat_Cov_DQ-022305.fasta", type = "DNA")
Deer_WT_OL <- read.FASTA("./Sequences/CoV_WTDeer_OL855841.fasta", type = "DNA")
Deer_WT_FJ <- read.FASTA("./Sequences/White-tailedDeerCov_FJ425187.1", type = "DNA")
Deer_Sambar <- read.FASTA("./Sequences/SambarDeerCov_FJ425189.1.fasta", type = "DNA")
Deer_Water <- read.FASTA("./Sequences/WaterDeer_MG518518.1.fasta", type = "DNA")

HCoV_229E <- read.FASTA("./Sequences/HCoV-229E_NC_002645", type = "DNA")
HCov_HK <- read.FASTA("./Sequences/HCoV-HKU1_NC_006577", type = "DNA")
HCov_NL <- read.FASTA("./Sequences/HCov-NL63_NC_005831", type = "DNA")
HCoV_OC <- read.FASTA("./Sequences/HCov-OC43_NC_006213", type = "DNA")
HCoV_BS <- read.FASTA("./Sequences/hCoV_BS001349.fasta", type = "DNA")
HCoV_ON <- read.FASTA("./Sequences/hCoV_ON078487.fasta", type = "DNA")
MERS <- read.FASTA("./Sequences/MERS-CoV_NC_019843", type = "DNA")
SARS_CoV <- read.FASTA("./Sequences/SARS-CoV_NC_004718", type = "DNA")
SARS_CoV2 <- read.FASTA("./Sequences/SARS-CoV2_NC_045512", type = "DNA")

# DNASbin to DNAStringSet
for(i in )

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


# Save Distance Matrix as .csv
write.csv(COVDMmat, "Cov_Distance.csv")
