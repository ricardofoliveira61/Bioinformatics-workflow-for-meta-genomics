# environment preparation
packages_to_check = c("BiocManager") 
installed_packages = installed.packages()[, "Package"]
missing_packages = setdiff(packages_to_check, installed_packages)

if (length(missing_packages) > 0) {
  message(paste("Installing missing packages:", paste(missing_packages, collapse = ", ")))
  install.packages(missing_packages, dependencies = TRUE)
} else {
  message("All required packages are already installed.")
}


library(BiocManager)

packages_to_check = c("DECIPHER","dada2","Biostrings") 
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


library(DECIPHER);message("DECHIPHER version: ",packageVersion("DECIPHER"))
library(dada2); message("DADA2 version: ",packageVersion("dada2"))
library(Biostrings); message("Biostrings version: ",packageVersion("Biostrings"))


# defining variables given by python
project = as.character(commandArgs(TRUE)[1])
workdir = as.character(commandArgs(TRUE)[2])
classifier = as.character(commandArgs(TRUE)[3])
fasta_file = as.character(commandArgs(TRUE)[4])
fasta_file = c("C:\\Users\\ricar\\Desktop\\Teste\\test\\Results\\test_asv_single_refseq.fa")
data_base = as.character(commandArgs(TRUE)[5])
load("~/tax/IDTaxa/SILVA_SSU_r132_March2018.RData") # CHANGE TO THE PATH OF YOUR TRAINING SET

setwd(workdir)
 
if(classifier == "DECIPHER"){
  # reading the ASVs fastafile
  dna = readDNAStringSet(fasta_file) 
  # using Decipher to classify each ASV
  ids = IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=TRUE)
  #ids_ = IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
  # ranks of interest
  ranks = c("domain", "phylum", "class", "order", "family", "genus", "species") 
  
  # Convert the output object of class "Taxa" to a matrix
  taxid = t(sapply(ids, function(x) {
    m = match(ranks, x$rank)
    taxa = x$taxon[m]
    taxa[startsWith(taxa, "unclassified_")] <- NA
    taxa
  }))
  colnames(taxid) = ranks; rownames(taxid) = dna@ranges@NAMES
  
  # extracting the confidence of classification for each ASV
  confidence = c()
  for(row in ids){
    confidence = c(confidence,round(min(row$confidence)))
  }
  # adding the co
  taxid= cbind(taxid,confidence)
  
} else{
  
}

# criar código para guardar os dados, com base na basededados utilizada
plot(ids, trainingSet)
