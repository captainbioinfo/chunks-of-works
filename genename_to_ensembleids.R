# Install Bioconductor and required packages if not already installed
install.packages("BiocManager")
BiocManager::install()

# Install biomaRt and clusterProfiler
BiocManager::install("biomaRt")
BiocManager::install("clusterProfiler")

# If an annotation package for Bos taurus is available, use it
# If not, you can try other approaches (e.g., GO annotations from Ensembl)

# Load necessary libraries
library(biomaRt)
library(clusterProfiler)

# Check if there's an available annotation database for Bos taurus
# org.Bt.eg.db for Bos taurus may not exist, so be sure to use the correct one
# Initialize mart
mart <- useMart("ENSEMBL_MART_ENSEMBL", host = "http://ensembl.org")

# List available datasets
datasets <- listDatasets(mart)

# Check the datasets available
head(datasets)

# Step 1: Connect to Ensembl Mart and list available datasets
mart <- useMart("ENSEMBL_MART_ENSEMBL")
datasets <- listDatasets(mart)

# Check the available datasets to find the correct one for Bos taurus
head(datasets)

# Correct dataset name for Bos taurus based on available datasets
# Example, if Bos taurus dataset is "bos_taurus_gene_ensembl"
mart <- useMart("ENSEMBL_MART_ENSEMBL", dataset = "btaurus_gene_ensembl")

# Step 2: Read gene names from CSV file
gene_list <- read.csv("k.csv")  # Replace with the actual path to your CSV file

# Check the structure of the data to ensure the correct column
head(gene_list)
gene_names <- gene_list$gene_name  # Adjust if necessary

# Retrieve Ensembl IDs for the gene names using biomaRt
gene_info <- getBM(attributes = c("external_gene_name", "ensembl_gene_id"),
                   filters = "external_gene_name", 
                   values = gene_names, 
                   mart = mart)

# Check the first few results to confirm the Ensembl IDs
head(gene_info)
write.csv(gene_info,"esemble_ids.csv")
################################################################################
