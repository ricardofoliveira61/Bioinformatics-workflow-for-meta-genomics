message("\nChecking if all the necessary R packages are installed")
packages_to_check = c("BiocManager", "ggplot2", "writexl", "openxlsx", "readxl","dplyr") 
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
packages_to_check = c("phyloseq") 
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


# loading packages
library(phyloseq); message("Phyloseq version: ",packageVersion("phyloseq"))
library(Biostrings)
library(dplyr)
library(ggplot2)


# creating R varivales from python inputs
project = as.character(commandArgs(TRUE)[1])
workdir = as.character(commandArgs(TRUE)[2])
asv_seqs = as.character(commandArgs(TRUE)[3])
asv_counts = as.character(commandArgs(TRUE)[4])
tax_class = as.character(commandArgs(TRUE)[5])
meta_samples = as.character(commandArgs(TRUE)[6])
setwd(workdir)
set.seed(123)


tax_class = "C:\\Users\\ricar\\Desktop\\Teste\\assignment_results\\teste_v4\\teste_RDP_classification_results.xlsx"
tax_class = "C:\\Users\\ricar\\Desktop\\Teste\\assignment_results\\teste_v4\\teste_SILVA_classification_results.xlsx"
asv_seqs = "C:\\Users\\ricar\\Desktop\\Teste\\V4\\Results\\V4_asv_single_refseq.fa"
asv_count = "C:\\Users\\ricar\\Desktop\\Teste\\V4\\Results\\V4_amplicon_results.xlsx"
meta_samples = "C:\\Users\\ricar\\Desktop\\Teste\\V4\\meta_mock1.csv"

read_data = function(filepath) {
  # Extract file extension
  extension = tools::file_ext(filepath)
  # Read data based on extension
  if (extension == "tsv") {
    read.table(filepath, header = TRUE, check.names = FALSE, sep = "\t")
  } else if (extension == "xlsx") {
    # Use readxl package to read xlsx files
    library(readxl)
    file = as.data.frame(read_excel(filepath, sheet = 1))
    rownames(file) = file[,1]; file = file[,-1]
  } else if(extension == "csv"){
    read.table(filepath, header = TRUE, row.names = 1, check.names = FALSE, sep = ",")
  } else    {
    stop("Unsupported file format. Please provide a tsv or xlsx file.")
  }
}


# tranformar os nomes das colunas em minusculas
#----------------------------------------------------------------------------------------------------------
# Data importation
# preparing the taxonomic classification dataframe
tax = as.data.frame(read_data(tax_class))
if ("confidence" %in% colnames(tax)){
  tax = tax[, !(colnames(tax) == "confidence")]
}
tax = as.matrix(tax)

# reading asv sequence file
seqs = Biostrings::readDNAStringSet(asv_seqs, format = "fasta")

# reading asv counts and preparing the dataframe
asv_counts = as.data.frame(read_data(asv_count))
counts = phyloseq::otu_table(asv_counts, taxa_are_rows = TRUE)

# reading metadata of the samples
meta = as.data.frame(read_data(meta_samples))


#----------------------------------------------------------------------------------------------------------
# creating phyloseq object
phylo = phyloseq(counts, seqs, sample_data(meta), tax_table(tax))
phylo


# extracting samples without taxonomic classification
phylo_no_domain = phylo %>% subset_taxa( domain == "NA" | is.na(domain))
phylo_no_domain


# extracting samples without phylum classification
phylo_no_phylum = phylo %>% subset_taxa( phylum == "NA" | is.na(phylum))
phylo_no_phylum


#EXport seqs for further inspection
#d_no_P %>%
#  refseq() %>%
#  Biostrings::writeXStringSet("./dada_not_class.fna", append=FALSE,
#                              compress=FALSE, compression_level=NA, format="fasta")


# remover os nas até à espéice
# removing the asvs withou classification
phylo_w_NA = phylo %>% subset_taxa (domain != "NA" | !is.na(domain))
phylo_w_NA


# selecting only the bacteria
phylo_bac_all = phylo_w_NA %>% subset_taxa( domain == "Bacteria" | is.na(domain))
phylo_bac_all


phylo_bac = filter_taxa(phylo_bac_all, function(x) {sum(x > 0) > 1}, TRUE)

# fazer com e sem normalização
# sem normalização (apenas filtragem)
#--------------------------------------------------------------------------------------------------------
# alpha diversity
plot_richness(phylo_bac , measures = c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher"),
              color = "kit" ) + ggtitle("Alpha diversity unnormalized")


#---------------------------------------------------------------------------------------------------------
# beta diversity
# atenção, os metadados precisas de ter pelo menos 2 variaveis, caso contrário bug
beta_div_PCoA = ordinate(phylo_bac, method = "PCoA", distance = "bray")
plot_ordination(phylo_bac, beta_div_PCoA, color = "kit") + ggtitle("Beta Diversity unnormalized PCoA")


beta_div_nmds = ordinate(phylo_bac, method="NMDS", distance = "bray")
plot_ordination(phylo_bac, beta_div_nmds, color = "kit") + ggtitle("Beta Diversity unnormalized NMDS")


#--------------------------------------------------------------------------------------------------------
# Normalization
# rarefaction
dbac_rare = rarefy_even_depth(phylo_bac, sample.size = min(sample_sums(phylo_bac)), replace = FALSE)

# Relative frequencies
dbacr = transform_sample_counts(phylo_bac, function(x) x/sum(x))


#--------------------------------------------------------------------------------------------------------
# alpha diversity normalized
plot_richness(dbac_rare , measures = c("Observed", "Chao1", "Shannon", "InvSimpson", "Fisher"),
              color = "kit") + ggtitle("Alpha diversity normalized")


#---------------------------------------------------------------------------------------------------------
# beta diversity mormalized
# atenção, os metadados precisas de ter pelo menos 2 variaveis, caso contrário bug
#PCoA
beta_div_PCoA_n = ordinate(dbac_rare, method = "PCoA", distance = "bray")
plot_ordination(dbac_rare, beta_div_PCoA_n, color = "kit") + ggtitle("Beta Diversity normalized PCoA")


#NMDS
beta_div_nmds_n = ordinate(dbac_rare, method="NMDS", distance = "bray")
plot_ordination(dbac_rare, beta_div_nmds_n, color = "kit") + ggtitle("Beta Diversity normalized NMDS")


#--------------------------------------------------------------------------------------------------------
#Taxonomic composition

phyloseq::plot_bar(dbac_rare, fill = "phylum") + ggtitle("Taxonomic composition: Phylum")

phyloseq::plot_bar(dbac_rare, fill = "class") + ggtitle("Taxonomic composition: Class")

phyloseq::plot_bar(dbac_rare, fill = "order") + ggtitle("Taxonomic composition: Order")

phyloseq::plot_bar(dbac_rare, fill = "family") + ggtitle("Taxonomic composition: Family")

phyloseq::plot_bar(dbac_rare, fill = "genus") + ggtitle("Taxonomic composition: Genus")


phylo_bac_tax_algo = phyloseq::tax_glom(dbac_rare, "phylum", NArm = TRUE)
phyloseq::plot_bar(phylo_bac_tax_algo, fill = "phylum") 

phylo_bac_tax_algo = phyloseq::tax_glom(dbac_rare, "class", NArm = TRUE)
phyloseq::plot_bar(phylo_bac_tax_algo, fill = "class")

phylo_bac_tax_algo = phyloseq::tax_glom(dbac_rare, "order", NArm = TRUE)
phyloseq::plot_bar(phylo_bac_tax_algo, fill = "order")

phylo_bac_tax_algo = phyloseq::tax_glom(dbac_rare, "family", NArm = TRUE)
phyloseq::plot_bar(phylo_bac_tax_algo, fill = "family")

phylo_bac_tax_algo = phyloseq::tax_glom(dbac_rare, "genus", NArm = TRUE)
phyloseq::plot_bar(phylo_bac_tax_algo, fill = "genus")


