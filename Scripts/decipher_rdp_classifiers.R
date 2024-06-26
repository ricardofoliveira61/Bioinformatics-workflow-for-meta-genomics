# environment preparation
# checking if the necessary R packages are installed
message("\nChecking if all the necessary R packages are installed")
packages_to_check = c("BiocManager", "writexl","openxlsx", "ggplot2","export") 
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


message("\n Loading the necessary packages.")
# loading the necessary packages
library(DECIPHER);message("DECHIPHER version: ",packageVersion("DECIPHER"))
library(dada2); message("DADA2 version: ",packageVersion("dada2"))
library(Biostrings); message("Biostrings version: ",packageVersion("Biostrings"))
library(writexl); library(openxlsx); library(ggplot2)


# defining variables given by python
project = as.character(commandArgs(TRUE)[1])
workdir = as.character(commandArgs(TRUE)[2])
classifier = as.character(commandArgs(TRUE)[3])
database = as.character(commandArgs(TRUE)[4])
fasta_file = as.character(commandArgs(TRUE)[5])
main_path = as.character(commandArgs(TRUE)[6])
parent_folder = dirname(main_path)
setwd(workdir)


message(paste0("\nLoading the database ",database))
# loading the training set of each database
if(database == "SILVA"){
  
  load(paste0(parent_folder, "/Databases/SILVA_SSU_r138_2019.RData"))

} else if (database == "RDP"){
  
  load(paste0(parent_folder, "/Databases/RDP_v18_July2020.RData"))

} else if (database == "GREENGENES"){
  
  print("Not avaible yet")
  q()
  
  }


# cheking if the directory assignment_results exist
assignment_path = paste0(workdir,"/assignment_results")
if (!file.exists(assignment_path)) {  
  dir.create(assignment_path)
}


# creating the project directory
assignment_path_results = paste0(assignment_path,"/",project)
if (!file.exists(assignment_path_results)) {  
  dir.create(assignment_path_results)
}


message(paste0("\nRunning the classifier ",classifier))
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
  taxid= cbind(taxid,confidence); taxid = as.data.frame(taxid)
  taxid$confidence = as.numeric(taxid$confidence)
  
  message("\nSaving the Taxonomic classification results")
  # Creating a workbook
  workbook = createWorkbook()
  
  # Adding DataFrames to separate sheets
  addWorksheet(workbook, sheetName = "Taxonomies")
  writeData(wb=workbook, sheet = workbook$sheetOrder[1], x = taxid, rowNames = T)
  
  # Saving the workbook
  saveWorkbook(workbook, file = paste0(assignment_path_results,"/", project,"_", database,
                                       "_classification_results.xlsx"), overwrite = TRUE)
  
  #plot(ids, trainingSet)
  #ggsave(filename = paste0(project,"_",database,"_classification_results.png"), 
  #     path = assignment_path_results)
  #width = 500, height = 500, units = "mm"
} else{
  
  print("Classifer not avaible yet")
  q()
  
  }
