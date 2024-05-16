library("BiocManager")
library("dada2"); packageVersion("dada2")
library(ggplot2)

# can be defined in python
setwd("C:/Users/ricar/Desktop/Git_hub/Projeto/Bioinformatics-workflow-for-meta-genomics")
# given on python
path_fastq = "./Data/filtered_mocks/Illumina/tests"

#creating the filterd_reads directory
filt_path = paste0(path_fastq,"/filtered_reads")
if (!file.exists(filt_path_test)) {
  dir.create(filt_path)}


# reading the forward and reverse reads files
illumina_fw = sort(list.files(path_fastq, pattern="_R1.fastq", full.names = TRUE))
illumina_rv = sort(list.files(path_fastq, pattern="_R2.fastq", full.names = TRUE))
# extracing samples names
sample_names = sapply(strsplit(basename(illumina_fw), "\\_"), `[`, 1)


# creating the directory to store the plots
plot_path= paste0(path_fastq,"/plots")
if (!file.exists(plot_path)) {  
  dir.create(plot_path)}


# creating the directory to store the quality profile plots of the unfiltered fastq files
plot_path_unfiltered= paste0(path_fastq,"/plots/quality_unfiltered")
if (!file.exists(plot_path_unfiltered)) {  
  dir.create(plot_path_unfiltered)}


for (I in 1:length(illumina_fw)) {
  #print(I)
  #print(illumina_fw[I])
  quality_plots = plotQualityProfile(illumina_fw[I])
  ggsave(filename = paste0(sample_names[I],"_qualityprofile_unfiltered_Fw.png"),
         path = plot_path_unfiltered)}


for (I in 1:length(illumina_rv)) {
  print(I)
  print(illumina_rv[I])
  quality_plots = plotQualityProfile(illumina_rv[I])
  ggsave(filename = paste0(sample_names[I],"_qualityprofile_unfiltered_Rv.png"),
         path = plot_path_unfiltered)}


# creating the path to store the filtered fastq files
filtFw = file.path(filt_path, paste0(sample_names, "_Fw_filt.fastq"))
filtRv = file.path(filt_path, paste0(sample_names, "_Rv_filt.fastq"))
names(filtFw) = sample_names
names(filtRv) = sample_names


# Filtration of fastq files
out = filterAndTrim(illumina_fw, filtFw, illumina_rv, filtRv, truncLen=c(240,240),
                    minLen=200, trimLeft=10, truncQ=2, maxEE=c(1,2), maxN=0, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE, multithread=TRUE)


illumina_fw_filt = sort(list.files(filt_path, pattern="_Fw_filt.fastq", full.names = TRUE))
illumina_rv_filt = sort(list.files(filt_path, pattern="_Rv_filt.fastq", full.names = TRUE))
# extracing sample names
sample_names_filt = sapply(strsplit(basename(illumina_fw_filt), "\\_"), `[`, 1)


#creating the directory to store the quality profile plots of the filtered fastq files
plot_path_filtered= paste0(path_fastq,"/plots/quality_filtered")
if (!file.exists(plot_path_filtered)) {  
  dir.create(plot_path_filtered)}


for (I in 1:length(illumina_fw_filt)) {
  #print(I)
  #print(illumina_fw_filt[I])
  quality_plots = plotQualityProfile(illumina_fw_filt[I])
  ggsave(filename = paste0(sample_names_filt[I],"_qualityprofile_filtered_Fw.png"),
         path = plot_path_filtered)}


for (I in 1:length(illumina_rv_filt)) {
  #print(I)
  #print(illumina_rv_filt[I])
  quality_plots = plotQualityProfile(illumina_rv_filt[I])
  ggsave(filename = paste0(sample_names_filt[I],"_qualityprofile_filtered_Rv.png"),
         path = plot_path_filtered)}


#estimate the error model for DADA2 algorithm using forward reads
errF = learnErrors(illumina_fw_filt, multithread=TRUE)
#estimate the error model for DADA2 algorithm using reverse reads
errR = learnErrors(illumina_rv_filt, multithread=TRUE)


plot_errors_path = paste0(path_fastq,"/plots/errors")
if (!file.exists(plot_errors_path)) {  
  dir.create(plot_errors_path)}


plotErrors(errF, nominalQ=TRUE)
ggsave(filename = "errors_Fw.png", path = plot_errors_path, width = 300, height = 300,
       units = "mm")


plotErrors(errR, nominalQ=TRUE)
ggsave(filename = "errors_Rv.png", path = plot_errors_path, width = 300, height = 300,
       units = "mm")


# ASV creation using the error models
dadaFs = dada(illumina_fw_filt, err=errF, multithread=TRUE) 

dadaRs = dada(illumina_rv_filt, err=errR, multithread=TRUE) 

merger = mergePairs(dadaFs, illumina_fw_filt, dadaRs,
                    illumina_rv_filt, verbose=TRUE) 


#table construction of all ASVs
seqtab = makeSequenceTable(merger)
# removing chimeras
seqtab_no_chim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#sum(seqtab_no_chim)/sum(seqtab) # perguntar o que isto faz

# filtering the results to only have the V4 amplicon Data
seqtab_V4 = seqtab_no_chim[,nchar(colnames(seqtab_no_chim)) %in% 230:240]


# removing singletons
is_1 = colSums(seqtab_V4) <= 1
seqtab_V4 = seqtab_V4[,!is_1]

#------------------------------------------------------------------------------------------------
# tracking the number of reads in each step of the workflow
# creating the function get N, returns the total number of unique sequences
getN = function(x) sum(getUniques(x))


# removing samples that are not used
#retriving the indices of the samples not used
idx_samples_not_used = which(out[,"reads.out"]==0)
if (length(idx_samples_not_used) !=0){
  # retriving the name of the samples not used
  # creating a matrix with the data of the samples not used
  samples_not_used = matrix(out[idx_samples_not_used,],nrow = length(idx_samples_not_used)
                            ,ncol = 2, byrow= TRUE)
  #giving to each row the name of each sample
  rownames(samples_not_used) = sample_names[idx_samples_not_used]
  # filtering the out matrix containing the results of read filtration
  out_filtered = out[-idx_samples_not_used,]
  # creating the missing columns
  new_columns = matrix(0, nrow = nrow(samples_not_used), ncol = 4)
  # adding the missing columns to the orinal matrix
  samples_not_used = cbind(samples_not_used, new_columns)
  # assigning the name of the columns 
  colnames(samples_not_used) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")

  } else {
  out_filtered = out
}


# creation of the track matrix
if(length(sample_names_filt) != 1){

  track = cbind(out_filtered, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(merger, getN), rowSums(seqtab_no_chim))

} else{
  
  track = cbind(out_filtered, getN(dadaFs), getN(dadaRs), getN(merger), rowSums(seqtab_no_chim))
}

colnames(track) = c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) = sample_names_filt

if(length(idx_samples_not_used) != 0){
  # adding the rows of the not used samples
  track = rbind(track,samples_not_used)
}


write.csv(track, "dada2_processing.csv")


#------------------------------------------------------------------------------------------------
# Data exportation
asv_dada = seqtab_V4

# giving our seq headers more manageable names (ASV1, ASV2...)
asv_seqs = colnames(asv_dada)
asv_headers = vector(dim(asv_dada)[2], mode="character")

for (i in 1:dim(asv_dada)[2]) {
  asv_headers[i] = paste(">ASV", i, sep="")
}

# making and writing out a fasta of our final ASV seqs:
asv_fasta = c(rbind(asv_headers, asv_seqs))
write(asv_fasta, "dada_single_refseq.fa")

# count table:
asv_tab = t(asv_dada)
row.names(asv_tab) = sub(">", "", asv_headers)
write.table(asv_tab, "dada_single_counts.tsv", sep="\t", quote=F, col.names=NA)
