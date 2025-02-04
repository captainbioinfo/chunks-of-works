install.packages("BiocManager")

BiocManager::install("DESeq2")

library(DESeq2)
#library(airway)
library(tidyverse)
# reads in count data
counts_data <- read.csv("gene_count_matrix.csv")
head(counts_data)
col_data <- read.csv("f4_timepoint2.csv")
all(colnames(counts_data) != rownames(col_data))

all(colnames(counts_data) == rownames(col_data))
all(sort(colnames(counts_data)) == sort(rownames(col_data)))
library(DESeq2)
#counts_data <- counts_data[, rownames(col_data)]
all(rownames(col_data) == colnames(counts_data))


# Remove the "gene_id" column from counts_data
counts_data <- counts_data[, !colnames(counts_data) %in% "gene_id"]

# Update the row names of col_data to match the column names of counts_data
rownames(col_data) <- colnames(counts_data)
 
# Check if column names of counts_data match row names of col_data
 all(sort(colnames(counts_data)) == sort(rownames(col_data)))  # Should return TRUE
 dds <- DESeqDataSetFromMatrix(countData =  counts_data,
                               colData = col_data,
                               design = ~ condition)
dds
# filtering low qulity reads
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

dds
dds$condition <- relevel(dds$condition, ref = "Control")

dds <- DESeq(dds)
res <- results(dds)
res
summary(res)
resultsNames(dds)
#results(dds, contrast = c("condition"," control","Treatment_24"))


# MA plot
plotMA(res)
############################
# Get all significant genes (padj < 0.05)
sig_genes <- subset(res, padj < 0.05)

# Get upregulated genes (log2FoldChange > 0 and padj < 0.05)
up_genes <- subset(res, padj < 0.05 & log2FoldChange > 0)

# Get downregulated genes (log2FoldChange < 0 and padj < 0.05)
down_genes <- subset(res, padj < 0.05 & log2FoldChange < 0)

# Convert to data frames and add gene names as a column
sig_genes_df <- as.data.frame(sig_genes)
sig_genes_df$gene_name <- rownames(sig_genes_df)

up_genes_df <- as.data.frame(up_genes)
up_genes_df$gene_name <- rownames(up_genes_df)

down_genes_df <- as.data.frame(down_genes)
down_genes_df$gene_name <- rownames(down_genes_df)

# Save to CSV files
write.csv(sig_genes_df, "significant_genes.csv", row.names = FALSE)
write.csv(up_genes_df, "upregulated_genes.csv", row.names = FALSE)
write.csv(down_genes_df, "downregulated_genes.csv", row.names = FALSE)

# Print summary of results
cat("Total significant genes:", nrow(sig_genes_df), "\n")
cat("Upregulated genes:", nrow(up_genes_df), "\n")
cat("Downregulated genes:", nrow(down_genes_df), "\n")

#######################################
# Install if you haven't already
#BiocManager::install("EnhancedVolcano")
#BiocManager::install("ggrepel")

# Load the package
library(EnhancedVolcano)
# For a more customized version with different colors and thresholds:
# First, let's make sure we have proper labels
# Convert results to a data frame and add gene names if they're not already present
res_df <- as.data.frame(res)
res_df$labels <- rownames(res_df)
EnhancedVolcano(res_df,
                lab = ifelse(res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, rownames(res_df), ""),
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Treatment vs Control',
                subtitle = 'Differential Expression Analysis',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4.0,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 0.4,
                legendPosition = 'right',
                legendLabSize = 10,
                drawConnectors = FALSE,
                widthConnectors = 0.5)


########

#####################
EnhancedVolcano(res_df,
                lab = "",  # Remove gene labels by setting this to an empty string
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Treatment vs Control',
                subtitle = 'Differential Expression Analysis',
                pCutoff = 0.05,
                FCcutoff = 1,
                pointSize = 3.0,
                labSize = 4.0,
                col = c("grey30", "forestgreen", "royalblue", "red2"),
                colAlpha = 0.4,
                legendPosition = 'right',
                legendLabSize = 10,
                drawConnectors = FALSE,
                widthConnectors = 0.5)


