#preparing R environment
message("\nChecking if all the necessary R packages are installed")
packages_to_check = c("BiocManager", "ggplot2", "writexl", "openxlsx") 
installed_packages = installed.packages()[, "Package"]
missing_packages = setdiff(packages_to_check, installed_packages)

if (length(missing_packages) > 0) {
  message(paste("Installing missing packages:", paste(missing_packages, collapse = ", ")))
  for(package in missing_packages){
    install.packages(package, dependencies = TRUE, repos = "https://cran.r-project.org")
  }
  
} else {
  message("All required packages are already installed.")
}


# checking if all BioConductor packages are installed
library(BiocManager)
packages_to_check = c("dada2") 
installed_packages = installed.packages()[, "Package"]
missing_packages = setdiff(packages_to_check, installed_packages)
if (length(missing_packages) > 0) {
  message(paste("Installing missing packages:", paste(missing_packages, collapse = ", ")))
  for (package in missing_packages){
    BiocManager::install(package)
  }
} else {
  message("All required BioConductor packages are already installed.")
}


#loading necessary libraries
library("dada2"); message("DADA2 version: ",packageVersion("dada2"))
library(ggplot2); library(writexl); library(openxlsx)


# creating R variables from python input
project = as.character(commandArgs(TRUE)[1])
workdir = as.character(commandArgs(TRUE)[2])
fastq_files = as.character(commandArgs(TRUE)[3])
region = as.character(commandArgs(TRUE)[4])
setwd(workdir)


# creating project folder
project_path = paste0(workdir,"/",project)
if (!file.exists(project_path)) {
  dir.create(project_path)
}


#creating the filterd_reads directory
filt_path = paste0(workdir,"/",project,"/filtered_reads")
if (!file.exists(filt_path)) {
  dir.create(filt_path)
}


# reading the forward and reverse reads files
illumina_fw = sort(list.files(fastq_files, pattern="_R1.fastq", full.names = TRUE))
illumina_rv = sort(list.files(fastq_files, pattern="_R2.fastq", full.names = TRUE))
# extracing samples names
sample_names = sapply(strsplit(basename(illumina_fw), "\\_"), `[`, 1)


# creating the directory to store the plots
plot_path= paste0(workdir,"/",project,"/plots")
if (!file.exists(plot_path)) {  
  dir.create(plot_path)
}


# creating the directory to store the quality profile plots of the unfiltered fastq files
plot_path_unfiltered= paste0(plot_path,"/quality_unfiltered")
if (!file.exists(plot_path_unfiltered)) {  
  dir.create(plot_path_unfiltered)
}


message("\nSaving the quality profile of the unfiltered forward reads")
for (I in 1:length(illumina_fw)) {
  #print(I)
  #print(illumina_fw[I])
  quality_plots = plotQualityProfile(illumina_fw[I])
  ggsave(filename = paste0(sample_names[I],"_qualityprofile_unfiltered_Fw.png"),
         path = plot_path_unfiltered)
}


message("\nSaving the quality profile of the unfiltered reverse reads")
for (I in 1:length(illumina_rv)) {
  #print(I)
  #print(illumina_rv[I])
  quality_plots = plotQualityProfile(illumina_rv[I])
  ggsave(filename = paste0(sample_names[I],"_qualityprofile_unfiltered_Rv.png"),
         path = plot_path_unfiltered)
}


# creating the path to store the filtered fastq files
filtFw = file.path(filt_path, paste0(sample_names, "_Fw_filt.fastq"))
filtRv = file.path(filt_path, paste0(sample_names, "_Rv_filt.fastq"))
names(filtFw) = sample_names
names(filtRv) = sample_names


message("\nRead filtering starting")
# Filtration of fastq files
if(region == "V4"){
out = filterAndTrim(illumina_fw, filtFw, illumina_rv, filtRv, truncLen=c(240,240),
                    minLen=200, trimLeft=10, truncQ=2, maxEE=c(1,2), maxN=0, rm.phix=TRUE,
                    compress=TRUE, verbose=FALSE, multithread=TRUE)
} else if(region=="V3-V4"){
  out = filterAndTrim(illumina_fw, filtFw, illumina_rv, filtRv, truncLen=c(249,249),
                    minLen=200, trimLeft=10, truncQ=2, maxEE=c(1,2), maxN=0, rm.phix=TRUE,
                    compress=TRUE, verbose=FALSE, multithread=TRUE)
}
message("\nRead filtering finished")


# extracting the path to the filtered fastq files
illumina_fw_filt = sort(list.files(filt_path, pattern="_Fw_filt.fastq", full.names = TRUE))
illumina_rv_filt = sort(list.files(filt_path, pattern="_Rv_filt.fastq", full.names = TRUE))


# extracing sample names
sample_names_filt = sapply(strsplit(basename(illumina_fw_filt), "\\_"), `[`, 1)



#creating the directory to store the quality profile plots of the filtered fastq files
plot_path_filtered= paste0(plot_path,"/quality_filtered")
if (!file.exists(plot_path_filtered)) {  
  dir.create(plot_path_filtered)
}


message("\nSaving the quality profile of the filtered forward reads")
for (I in 1:length(illumina_fw_filt)) {
  #print(I)
  #print(illumina_fw_filt[I])
  quality_plots = plotQualityProfile(illumina_fw_filt[I])
  ggsave(filename = paste0(sample_names_filt[I],"_qualityprofile_filtered_Fw.png"),
         path = plot_path_filtered)
}


message("\nSaving the quality profile of the filtered reverse reads")
for (I in 1:length(illumina_rv_filt)) {
  #print(I)
  #print(illumina_rv_filt[I])
  quality_plots = plotQualityProfile(illumina_rv_filt[I])
  ggsave(filename = paste0(sample_names_filt[I],"_qualityprofile_filtered_Rv.png"),
         path = plot_path_filtered)
}


message("\nLearning the error model for the forward and reverse reads")
#estimate the error model for DADA2 algorithm using forward reads
errF = learnErrors(illumina_fw_filt, multithread=TRUE)
#estimate the error model for DADA2 algorithm using reverse reads
errR = learnErrors(illumina_rv_filt, multithread=TRUE)


plot_errors_path = paste0(plot_path,"/errors")
if (!file.exists(plot_errors_path)) {  
  dir.create(plot_errors_path)}


message("\nSaving the errors plots")
plotErrors(errF, nominalQ=TRUE)
ggsave(filename = "errors_Fw.png", path = plot_errors_path, width = 300, height = 300,
       units = "mm")


plotErrors(errR, nominalQ=TRUE)
ggsave(filename = "errors_Rv.png", path = plot_errors_path, width = 300, height = 300,
       units = "mm")


message("\nCreating the ASVs")
# ASV creation using the error models
dadaFs = dada(illumina_fw_filt, err=errF, multithread=TRUE) 

dadaRs = dada(illumina_rv_filt, err=errR, multithread=TRUE) 

message("\nMerging the forward and reverse ASVs")
merger = mergePairs(dadaFs, illumina_fw_filt, dadaRs,
                    illumina_rv_filt, verbose=TRUE) 


message("\nProcessing ASVs results")
#table construction of all ASVs
seqtab = makeSequenceTable(merger)
#table(nchar(getSequences(seqtab))) 
# removing chimeras
seqtab_no_chim = removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
#table(nchar(getSequences(seqtab_no_chim))) 
#sum(seqtab_no_chim)/sum(seqtab) # perguntar o que isto faz

# filtering the results to only have the V4 amplicon Data or V3-V4
if(region=="V4"){
  seqtab_filt = seqtab_no_chim[,nchar(colnames(seqtab_no_chim)) %in% 230:240]  
} else if(region== "V3-V4") {
  seqtab_filt = seqtab_no_chim[,nchar(colnames(seqtab_no_chim)) %in% 420:450]
}


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


#------------------------------------------------------------------------------------------------
# Data exportation
# creating the results directory
results_path = paste0(workdir,"/",project,"/Results")
if (!file.exists(results_path)) {
  dir.create(results_path)}


asvs_dada = seqtab_filt

# giving our seq headers more manageable names (ASV1, ASV2...)
asv_seqs = colnames(asvs_dada)
names_asv = vector(dim(asvs_dada)[2], mode="character")

for (i in 1:dim(asvs_dada)[2]) {
  names_asv[i] = paste(">ASV", i, sep="")
}


message("\nSaving the results in a fasta file")
# making and writing out a fasta of the final ASV seqs:
asv_fasta = c(rbind(names_asv, asv_seqs))
write(asv_fasta, paste0(results_path,"/",project,"_asv_single_refseq.fa"))

# count table:
asv_tab = as.data.frame(t(asvs_dada))
row.names(asv_tab) = sub(">", "", names_asv)
colnames(asv_tab) = sample_names_filt


message("\nSaving analysis status and ASVs count into a excel file")
# Creating a workbook
workbook = createWorkbook()

# Adding DataFrames to separate sheets
addWorksheet(workbook, sheetName = "Number of reads each step")
writeData(wb=workbook, sheet = workbook$sheetOrder[1], x = track, rowNames = T)

addWorksheet(workbook, sheetName = "ASVs count")
writeData(wb=workbook, sheet = workbook$sheetOrder[2], x = asv_tab, rowNames = T)

# Saving the workbook
saveWorkbook(workbook, file = paste0(results_path,"/", project,"_amplicon_results.xlsx"), overwrite = TRUE)
