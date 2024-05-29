# environment preparation
# cheking if the necessary R packages are installed

message("\nCheking if all the necessary R packages are installed")
packages_to_check = c("BiocManager") 
installed_packages = installed.packages()[, "Package"]
missing_packages = setdiff(packages_to_check, installed_packages)
if (length(missing_packages) > 0) {
  message(paste("Installing missing packages:", paste(missing_packages, collapse = ", ")))
  install.packages(missing_packages, dependencies = TRUE)
} else {
  message("All required packages are already installed.")
}


# checking if the necessary BioConductor packages are installed
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


# loading the necessary packages
library(DECIPHER);message("DECHIPHER version: ",packageVersion("DECIPHER"))
library(dada2); message("DADA2 version: ",packageVersion("dada2"))
library(Biostrings); message("Biostrings version: ",packageVersion("Biostrings"))


# defining variables given by python
project = as.character(commandArgs(TRUE)[1])
workdir = as.character(commandArgs(TRUE)[2])
classifier = as.character(commandArgs(TRUE)[3])
database = as.character(commandArgs(TRUE)[4])
fasta_file = as.character(commandArgs(TRUE)[5])
main_path = as.character(commandArgs(TRUE)[6])
parent_folder = dirname(main_path)


# changing the working directory
setwd(workdir)

message(paste0("\nLoading the database",database))
# loading the training set of each database
if(database == "SILVA"){
  load(paste0(parent_folder, "/Databases/SILVA_SSU_r138_2019.RData"))
  print("SILVA")
} else if (database == "RDP"){
  load(paste0(parent_folder, "/Databases/RDP_v18_July2020.RData"))
  print("RDP")
} else if (database == "GREENGENES"){
  print("GREEN")
  }


# cheking if the directory assignment_results exist
assignment_path = paste0(workdir,"/assignment_results")
if (!file.exists(assignment_path)) {  
  dir.create(assignment_path)
  }


message(paste0("\nRunning the classifier",classifier))
# running the classifiers
if(classifier == "DECIPHER"){

  dna = readDNAStringSet(fasta_file) 
  ids = IdTaxa(dna, trainingSet, strand="both", processors=NULL, verbose=TRUE)
  #ids_ = IdTaxa(dna, trainingSet, strand="top", processors=NULL, verbose=FALSE)
  ranks = c("domain", "phylum", "class", "order", "family", "genus", "species") 
  
  # Converting the output to a more readable form
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
  
  # adding the cofindece to the matrix
  taxid= cbind(taxid,confidence)
  
} else{
  print("Classifer not avaible yet")
  }

# criar código para guardar os dados, com base na basededados utilizada
plot(ids, trainingSet)


