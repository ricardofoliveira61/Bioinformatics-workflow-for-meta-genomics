# preparing the environment 
setwd("C:/Users/ricar/Desktop/Git_hub/Projeto/Bioinformatics-workflow-for-meta-genomics")
getwd()

# installation of DADA2
#BiocManager::install("dada2",force = T)
#

# loading of the required packages
library("BiocManager")
library("dada2"); packageVersion("dada2")
library(ggplot2)


# ONT amplicon data
path = "./Data/filtered_mocks/ONT"
list.files(path)
ad_ont = sort(list.files(path, pattern="DRR223274.fastq", full.names = TRUE))
sample_names = sapply(strsplit(basename(ad_ont), "\\."), `[`, 1)

plotQualityProfile(ad_ont[1])



# Illumina amplicon data
# for a selected set of test mocks

path_test = "./Data/filtered_mocks/Illumina/tests"
#creating the filterd_reads directory
filt_path_test = paste0(path_test,"/filtered_reads")
if (!file.exists(filt_path_test)) {
  dir.create(filt_path_test)}

# reading the forward and reverse read files for Zymmo and ATCT mocks
ad_illumina_fw = sort(list.files(path_test, pattern="_R1.fastq", full.names = TRUE))
ad_illumina_rv = sort(list.files(path_test, pattern="_R2.fastq", full.names = TRUE))


# extracing sample names
sample_names_test = sapply(strsplit(basename(ad_illumina_fw), "\\_"), `[`, 1)


# ploting the quality of the first 3 foward sammples
if (!file.exists("./Data/filtered_mocks/Illumina/plots")) {  # Check if directory exists (optional)
  dir.create("./Data/filtered_mocks/Illumina/plots")  # Create directory if needed (optional)
}
plotQualityProfile(ad_illumina_fw[1:3])
ggplot() # adding the title
title("Quality profile")
# saving the plot
ggsave(filename = "quality_unfiltered.png", path = paste0(path_test,"/plots"),
       width = 496 , height = 443, limitsize = F, units = "mm")


# creating the path to store the filtered fastq files
filtFw = file.path(filt_path_test, paste0(sample_names_test, "_Fw_filt.fastq"))
filtRv = file.path(filt_path_test, paste0(sample_names_test, "_Rv_filt.fastq"))
names(filtFw) = sample_names_test
names(filtRv) = sample_names_test


# using the funcion to filter and trim the fastQ files 
# using windows can take longer to filter since it the pacakge can't use multithread
out = filterAndTrim(ad_illumina_fw, filtFw, ad_illumina_rv, filtRv, truncLen=c(240,240),
                    minLen=200, trimLeft=10, truncQ=2, maxEE=c(1,2), maxN=0, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE, multithread=TRUE)

ad_illumina_fw_filt = sort(list.files(filt_path_test, pattern="_Fw_filt.fastq", full.names = TRUE))
ad_illumina_rv_filt = sort(list.files(filt_path_test, pattern="_Rv_filt.fastq", full.names = TRUE))
# extracing sample names
sample_names_filt = sapply(strsplit(basename(ad_illumina_fw_filt), "\\_"), `[`, 1)

# ploting the quality profile of the filtered fastQ files
plotQualityProfile(filtFw[1:3])
title("quality profile")
# saving the plot
ggsave(filename = "quality_filtered.png", path = paste0(path_test,"/plots"),
       width = 496 , height = 443, limitsize = F, units = "mm")


#estimate the error model for DADA2 algorithm using forward reads
errF = learnErrors(ad_illumina_fw_filt, multithread=TRUE)
#estimate the error model for DADA2 algorithm using reverse reads
errR = learnErrors(ad_illumina_rv_filt, multithread=TRUE)


plotErrors(errF, nominalQ=TRUE)
ggsave(filename = "errors_fw.png", path = paste0(path_test,"/plots"),
       width = 496 , height = 443, limitsize = F, units = "mm")

plotErrors(errR, nominalQ=TRUE)
ggsave(filename = "errors_rv.png", path = paste0(path_test,"/plots"),
       width = 496 , height = 443, limitsize = F, units = "mm")


# ASV creation using the error models
dadaFs = dada(ad_illumina_fw_filt, err=errF, multithread=TRUE) #filtered Fw

dadaRs = dada(ad_illumina_rv_filt, err=errR, multithread=TRUE) # filtered Rv

merger = mergePairs(dadaFs, ad_illumina_fw_filt, dadaRs,
                    ad_illumina_rv_filt, verbose=TRUE) # merging fw and rv


#table construction of all ASVs
seqtab = makeSequenceTable(merger)
dim(seqtab)


# checking the number of ASVs with certain length 
table(nchar(getSequences(seqtab)))


# removing chimeras
seqtab_no_chim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_no_chim)
table(nchar(getSequences(seqtab_no_chim)))
sum(seqtab_no_chim)/sum(seqtab)

seqtab2 = seqtab_no_chim[,nchar(colnames(seqtab_no_chim)) %in% 231:240]
table(nchar(getSequences(seqtab2)))


# tracking the number of reads in each step of the workflow
getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merger, getN), rowSums(seqtab_no_chim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample_names_test
head(track)

# add the code to export the data

# for all Illumina samples
path = "./Data/filtered_mocks/Illumina"
#creating the filterd_reads directory
filt_path = "./Data/filtered_mocks/Illumina/filtered_reads"
if (!file.exists(filt_path)) {
  dir.create(filt_path)}

# reading the forward and reverse read files for Zymmo and ATCT mocks
ad_illumina_fw = sort(list.files(path, pattern="_1.fastq", full.names = TRUE))
ad_illumina_rv = sort(list.files(path, pattern="_2.fastq", full.names = TRUE))


# extracing sample names
sample_names = sapply(strsplit(basename(ad_illumina_fw), "\\_"), `[`, 1)


# ploting the quality of the first 3 foward sammples
if (!file.exists("./Data/filtered_mocks/Illumina/plots")) {  # Check if directory exists (optional)
  dir.create("./Data/filtered_mocks/Illumina/plots")  # Create directory if needed (optional)
}
plotQualityProfile(ad_illumina_fw[1:3])
# saving the plot
ggsave(filename = "quality_unfiltered.png", path = "./Data/filtered_mocks/Illumina/plots",
       width = 496 , height = 443, limitsize = F, units = "mm")




# creating the path to store the filtered fastq files
filtFw = file.path(filt_path, paste0(sample_names, "_Fw_filt.fastq"))
filtRv = file.path(filt_path, paste0(sample_names, "_Rv_filt.fastq"))
names(filtFw) = sample_names
names(filtRv) = sample_names


# using the funcion to filter and trim the fastQ files 
# using windows can take longer to filter since it the pacakge can't use multithread
out = filterAndTrim(ad_illumina_fw, filtFw, ad_illumina_rv, filtRv, truncLen=c(240,240),
                    minLen=200, trimLeft=10, truncQ=2, maxEE=2, maxN=0, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE, multithread=TRUE)

# ploting the quality profile of the filtered fastQ files
plotQualityProfile(filtFw[1:3])
# saving the plot
ggsave(filename = "quality_filtered.png", path = "./Data/filtered_mocks/Illumina/plots",
       width = 496 , height = 443, limitsize = F, units = "mm")


ad_illumina_fw_filt = sort(list.files(filt_path, pattern="_Fw_filt.fastq", full.names = TRUE))
ad_illumina_rv_filt = sort(list.files(filt_path, pattern="_Rv_filt.fastq", full.names = TRUE))


# extracing sample names
sample_names_filt = sapply(strsplit(basename(ad_illumina_fw), "\\_"), `[`, 1)



#estimate the error model for DADA2 algorithm using forward reads
errF = learnErrors(ad_illumina_fw_filt, multithread=TRUE)
#estimate the error model for DADA2 algorithm using reverse reads
errR = learnErrors(ad_illumina_rv_filt, multithread=TRUE)


plotErrors(errF, nominalQ=TRUE)
ggsave(filename = "errors_fw.png", path = "./Data/filtered_mocks/Illumina/plots",
       width = 496 , height = 443, limitsize = F, units = "mm")

plotErrors(errR, nominalQ=TRUE)
ggsave(filename = "errors_rv.png", path = "./Data/filtered_mocks/Illumina/plots",
       width = 496 , height = 443, limitsize = F, units = "mm")


# ASV creation using the error models
dadaFs = dada(ad_illumina_fw_filt, err=errF, multithread=TRUE) #filtered Fw

dadaRs = dada(ad_illumina_rv_filt, err=errR, multithread=TRUE) # filtered Rv

merger = mergePairs(dadaFs, ad_illumina_fw_filt, dadaRs,
                    ad_illumina_rv_filt, verbose=TRUE) # merging fw and rv


#table construction of all ASVs
seqtab = makeSequenceTable(merger)
dim(seqtab)


# checking the number of ASVs with certain length 
table(nchar(getSequences(seqtab)))


# removing chimeras
seqtab_no_chim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab_no_chim)
table(nchar(getSequences(seqtab_no_chim)))
sum(seqtab_no_chim)/sum(seqtab)



seqtab2 = seqtab[,nchar(colnames(seqtab)) %in% 250:256]

out
# tracking the number of reads in each step of the workflow
getN = function(x) sum(getUniques(x))
track = cbind(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merger, getN), rowSums(seqtab_no_chim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample_names_filt
head(track)


## pipline professor
# os datasets tem de ser varridos previamente para verificar se 
# primers estão presentes ou não nas extremidades 5' 
# utilizar pro exemplo o cutadpt

# para descarregar os SRA utiliza-se a ferramenta SRA tools
# sra formato nativo de armazenamento de dados de sequenciação
#usado pelo ncbi
# criar uma rotina para de forma automática analisar a presença dos primers



