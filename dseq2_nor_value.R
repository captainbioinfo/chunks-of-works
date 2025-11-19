# ==============================================================================
# DESeq2 Differential Gene Expression Analysis: Treatment vs Control
# ==============================================================================
install.packages("DESEq2")
install.packages("tidyverse")
# Set seed for reproducibility
set.seed(123)

# ==============================================================================
# 1. LOAD REQUIRED LIBRARIES
# ==============================================================================
cat("Loading required libraries...\n")
library(DESeq2)
library(tidyverse)

# ==============================================================================
# 2. LOAD DATA
# ==============================================================================
cat("Loading count data and sample metadata...\n")

counts_data <- read.csv("PRJNA551141_edited.CSV")
col_data <- read.csv("sampleinfo141.csv")

# Fix column names
gene_ids <- counts_data$GeneId
counts_data_no_gene_id <- counts_data[, -1]
colnames(col_data) <- c("SampleName", "condition")

# Match sample order
if (!all(colnames(counts_data_no_gene_id) == col_data$SampleName)) {
  cat("Reordering sample info to match count matrix...\n")
  col_data <- col_data[match(colnames(counts_data_no_gene_id), col_data$SampleName), ]
}

# ==============================================================================
# 3. CREATE DESeq2 OBJECT
# ==============================================================================
cat("Creating DESeq2 dataset...\n")
rownames(col_data) <- col_data$SampleName

dds <- DESeqDataSetFromMatrix(countData = counts_data_no_gene_id,
                              colData = col_data,
                              design = ~ condition)

# ==============================================================================
# 4. FILTER LOW COUNT GENES
# ==============================================================================
cat("Filtering low count genes...\n")
keep <- rowSums(counts(dds)) >= 0  # optional threshold
dds_filtered <- dds[keep, ]
gene_ids_filtered <- gene_ids[keep]

cat(sprintf("Filtered from %d to %d genes\n", nrow(dds), nrow(dds_filtered)))

# ==============================================================================
# 5. RUN DE ANALYSIS
# ==============================================================================
cat("Running DESeq2 analysis...\n")
dds_filtered <- DESeq(dds)

# ==============================================================================
# 6. GET RESULTS
# ==============================================================================
cat("Extracting results...\n")
res <- results(dds_filtered, contrast = c("condition", "Treatment", "Control"))
res_df <- as.data.frame(res)
res_df$gene_id <- gene_ids_filtered

# Get normalized counts
norm_counts <- counts(dds_filtered, normalized = TRUE)
norm_counts_df <- as.data.frame(norm_counts)
norm_counts_df$gene_id <- gene_ids_filtered

# ==============================================================================
# 7. FILTER SIGNIFICANT GENES
# ==============================================================================
padj_threshold <- 0.05
log2fc_threshold <- 0

cat("Filtering significant genes...\n")
up_genes <- subset(res_df, padj < padj_threshold & log2FoldChange > log2fc_threshold & !is.na(padj))
down_genes <- subset(res_df, padj < padj_threshold & log2FoldChange < -log2fc_threshold & !is.na(padj))

cat(sprintf("Upregulated: %d genes\n", nrow(up_genes)))
cat(sprintf("Downregulated: %d genes\n", nrow(down_genes)))

# ==============================================================================
# 8. FORMAT OUTPUT
# ==============================================================================
format_output_df_with_norm <- function(df, description) {
  cat(sprintf("Formatting %s (%d genes)...\n", description, nrow(df)))
  
  if (nrow(df) == 0) {
    return(data.frame())
  }
  
  norm_subset <- norm_counts_df %>%
    filter(gene_id %in% df$gene_id)
  
  df_formatted <- df %>%
    select(gene_id, baseMean, log2FoldChange, lfcSE, stat, pvalue, padj) %>%
    left_join(norm_subset, by = "gene_id") %>%
    arrange(padj) %>%
    mutate(across(where(is.numeric), ~ round(.x, 6)))
  
  return(df_formatted)
}

up_genes_formatted <- format_output_df_with_norm(up_genes, "upregulated genes")
down_genes_formatted <- format_output_df_with_norm(down_genes, "downregulated genes")

# ==============================================================================
# 9. SAVE RESULTS
# ==============================================================================
write.csv(up_genes_formatted, "upregulated_genes.csv", row.names = FALSE)
write.csv(down_genes_formatted, "downregulated_genes.csv", row.names = FALSE)

cat("Files saved:\n")
cat("- upregulated_genes.csv\n")
cat("- downregulated_genes.csv\n")
cat("✓ DE analysis completed successfully\n")

