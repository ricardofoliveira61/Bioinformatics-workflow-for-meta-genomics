# preparing the environment 
setwd("C:/Users/ricar/Desktop/Git_hub/Projeto/Bioinformatics-workflow-for-meta-genomics")

# installation of DADA2
BiocManager::install("dada2",force = T)
BiocManager::install("stats",force = T)


# loading of the required packages
library("BiocManager")
library("dada2"); packageVersion("dada2")


# ONT amplicon data
path = "./Data/filtered_mocks/ONT"
list.files(path)
ad_ont = sort(list.files(path, pattern="DRR223274.fastq", full.names = TRUE))
sample_names = sapply(strsplit(basename(ad_ont), "\\."), `[`, 1)

plotQualityProfile(ad_ont[1])
# the graphic was empty????

# Illumina amplicon data
path = "./Data/filtered_mocks/Illumina"
list.files(path)
ad_illumina = sort(list.files(path, pattern="SRR18716709.fastq", full.names = TRUE))
sample_names = sapply(strsplit(basename(ad_illumina), "\\."), `[`, 1)

plotQualityProfile(ad_illumina[2])
## this graph showed an back line???

