# ==============================================================================
# edgeR Differential Gene Expression Analysis
# Analysis: Treatment vs Control
# ==============================================================================

# Set seed for reproducibility
set.seed(123)

# ==============================================================================
# 1. LOAD REQUIRED LIBRARIES
# ==============================================================================
cat("Loading required libraries...\n")
library(edgeR)
library(tidyverse)
library(limma)
library(RColorBrewer)
library(pheatmap)

# For volcano plots (install if not available)
if (!require("EnhancedVolcano", quietly = TRUE)) {
  cat("Installing EnhancedVolcano for volcano plots...\n")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("EnhancedVolcano")
  library(EnhancedVolcano)
} else {
  library(EnhancedVolcano)
}

# ==============================================================================
# 2. DATA LOADING AND PREPARATION
# ==============================================================================
cat("Loading count data and sample metadata...\n")

# Load count matrix and sample information
counts_data <- read.csv("gene_count_matrix.csv")
col_data <- read.csv("pro_42_sampleinfo.csv")

# Store gene IDs before removing them from counts data
gene_ids <- counts_data$gene_id

# Remove gene_id column for edgeR processing (edgeR expects numeric matrix)
counts_data_no_gene_id <- counts_data[, -1]

# Convert to matrix (required for edgeR)
counts_matrix <- as.matrix(counts_data_no_gene_id)

# ==============================================================================
# 3. SAMPLE METADATA PREPARATION
# ==============================================================================
cat("Preparing sample metadata...\n")

# Display sample names for comparison
cat("Sample names in count data (first 10):\n")
print(head(colnames(counts_matrix), 10))
cat("\nSample names in metadata (first 10):\n")
print(head(rownames(col_data), 10))

# Check if samples match
samples_match <- all(colnames(counts_matrix) == rownames(col_data$SampleName))

if (!samples_match) {
  cat("\n❌ Sample names do not match! Attempting to fix...\n")
  
  # Check if col_data has a sample name column instead of rownames
  cat("Columns in metadata:\n")
  print(colnames(col_data))
  
  # Common fixes:
  # 1. Check if there's a sample ID column in col_data
  if ("SampleName" %in% colnames(col_data)) {
    cat("Found 'SampleName' column. Setting as rownames...\n")
    rownames(col_data) <- col_data$SampleName
  } else if ("sample_id" %in% colnames(col_data)) {
    cat("Found 'sample_id' column. Setting as rownames...\n")
    rownames(col_data) <- col_data$sample_id
  } else if ("Sample" %in% colnames(col_data)) {
    cat("Found 'Sample' column. Setting as rownames...\n")
    rownames(col_data) <- col_data$Sample
  } else if ("SampleID" %in% colnames(col_data)) {
    cat("Found 'SampleID' column. Setting as rownames...\n")
    rownames(col_data) <- col_data$SampleID
  } else {
    # If first column looks like sample names, use it
    cat("Using first column as sample names...\n")
    rownames(col_data) <- col_data[,1]
  }
  
  # 2. Check if samples are in different order
  count_samples <- colnames(counts_matrix)
  meta_samples <- rownames(col_data)
  
  cat(sprintf("Count data has %d samples\n", length(count_samples)))
  cat(sprintf("Metadata has %d samples\n", length(meta_samples)))
  
  # Find common samples
  common_samples <- intersect(count_samples, meta_samples)
  cat(sprintf("Common samples: %d\n", length(common_samples)))
  
  if (length(common_samples) == 0) {
    cat("❌ No common samples found!\n")
    cat("Count data samples:\n")
    print(count_samples)
    cat("Metadata samples:\n")
    print(meta_samples)
    stop("Cannot proceed - no matching samples between count data and metadata!")
  }
  
  if (length(common_samples) < length(count_samples)) {
    cat("⚠️  Warning: Some samples are missing from metadata\n")
    missing_in_meta <- setdiff(count_samples, meta_samples)
    cat("Missing in metadata:\n")
    print(missing_in_meta)
  }
  
  if (length(common_samples) < length(meta_samples)) {
    cat("⚠️  Warning: Some samples are missing from count data\n")
    missing_in_counts <- setdiff(meta_samples, count_samples)
    cat("Missing in count data:\n")
    print(missing_in_counts)
  }
  
  # Subset both datasets to common samples and ensure same order
  cat("Subsetting to common samples and reordering...\n")
  counts_matrix <- counts_matrix[, common_samples]
  col_data <- col_data[common_samples, , drop = FALSE]
  
  # Final check
  if (all(colnames(counts_matrix) == rownames(col_data))) {
    cat("✅ Sample order is now consistent!\n")
  } else {
    cat("❌ Still not matching after attempts to fix\n")
    stop("Cannot resolve sample name mismatch!")
  }
} else {
  cat("✅ Sample order is consistent!\n")
}

# ==============================================================================
# 4. CREATE edgeR DGEList OBJECT
# ==============================================================================
cat("Creating edgeR DGEList object...\n")

# Create DGEList object
dge <- DGEList(counts = counts_matrix, 
               samples = col_data,
               genes = data.frame(gene_id = gene_ids))

# Display basic information
cat(sprintf("DGEList created with %d genes and %d samples\n", 
            nrow(dge), ncol(dge)))
cat("Sample information:\n")
print(table(dge$samples$condition))

# ==============================================================================
# 5. QUALITY CONTROL AND FILTERING
# ==============================================================================
cat("Applying quality control filters...\n")

# Calculate library sizes and normalization factors
dge$samples$lib.size <- colSums(dge$counts)
cat("Library sizes:\n")
print(dge$samples$lib.size)

# Filter low-count genes
# Keep genes with at least 1 CPM in at least 2 samples (minimum group size)
min_samples <- min(table(col_data$condition))
cat(sprintf("Minimum group size: %d samples\n", min_samples))

keep <- filterByExpr(dge, group = col_data$condition, min.count = 10)
dge_filtered <- dge[keep, , keep.lib.sizes = FALSE]

cat(sprintf("Filtered from %d to %d genes\n", nrow(dge), nrow(dge_filtered)))

# Recalculate library sizes after filtering
dge_filtered$samples$lib.size <- colSums(dge_filtered$counts)

# ==============================================================================
# 6. NORMALIZATION
# ==============================================================================
cat("Calculating normalization factors...\n")

# Calculate normalization factors (TMM - Trimmed Mean of M-values)
dge_filtered <- calcNormFactors(dge_filtered, method = "TMM")

cat("Normalization factors:\n")
print(dge_filtered$samples$norm.factors)

# ==============================================================================
# 7. DATA EXPLORATION AND QUALITY ASSESSMENT
# ==============================================================================
cat("Performing data exploration...\n")

# Calculate log2 CPM for visualization
lcpm <- cpm(dge_filtered, log = TRUE)

# Create MDS plot (Multi-Dimensional Scaling)
png("MDS_plot.png", width = 800, height = 600, res = 300)
plotMDS(dge_filtered, 
        col = as.numeric(as.factor(col_data$condition)), 
        pch = 19, cex = 1.5,
        main = "MDS Plot: Sample Clustering")
legend("topright", 
       levels(as.factor(col_data$condition)), 
       col = 1:length(levels(as.factor(col_data$condition))), 
       pch = 19, cex = 0.8)
dev.off()
cat("MDS plot saved as 'MDS_plot.png' ✓\n")

# Create boxplot of log2 CPM
png("Boxplot_logCPM.png", width = 1000, height = 600, res = 300)
boxplot(lcpm, 
        las = 2, 
        main = "Log2 CPM Distribution Across Samples",
        ylab = "Log2 CPM",
        col = rainbow(ncol(lcpm)))
dev.off()
cat("Boxplot saved as 'Boxplot_logCPM.png' ✓\n")

# ==============================================================================
# 8. EXPERIMENTAL DESIGN MATRIX
# ==============================================================================
cat("Setting up experimental design...\n")

# Create design matrix
design <- model.matrix(~ 0 + condition, data = col_data)
colnames(design) <- levels(as.factor(col_data$condition))

cat("Design matrix:\n")
print(design)

# Define contrast matrix (Treatment vs Control)
contrast_matrix <- makeContrasts(
  TreatmentVsControl = Treatment - Control,
  levels = design
)

cat("Contrast matrix:\n")
print(contrast_matrix)

# ==============================================================================
# 9. ESTIMATE DISPERSIONS
# ==============================================================================
cat("Estimating dispersions...\n")

# Estimate common and tagwise dispersions
dge_filtered <- estimateDisp(dge_filtered, design)

cat(sprintf("Common dispersion: %.4f\n", dge_filtered$common.dispersion))
cat(sprintf("Trended dispersion range: %.4f - %.4f\n", 
            min(dge_filtered$trended.dispersion), 
            max(dge_filtered$trended.dispersion)))

# Plot biological coefficient of variation
png("BCV_plot.png", width = 800, height = 600, res = 300)
plotBCV(dge_filtered, main = "Biological Coefficient of Variation")
dev.off()
cat("BCV plot saved as 'BCV_plot.png' ✓\n")

# ==============================================================================
# 10. DIFFERENTIAL EXPRESSION TESTING
# ==============================================================================
cat("Performing differential expression testing...\n")

# Fit generalized linear model
fit <- glmQLFit(dge_filtered, design)

# Perform quasi-likelihood F-test
qlf <- glmQLFTest(fit, contrast = contrast_matrix)

# Extract results
results <- topTags(qlf, n = Inf, sort.by = "PValue")
results_df <- as.data.frame(results$table)

cat(sprintf("Total genes tested: %d\n", nrow(results_df)))

# ==============================================================================
# 11. DEFINE SIGNIFICANCE THRESHOLDS AND IDENTIFY SIGNIFICANT GENES
# ==============================================================================
cat("Setting significance thresholds...\n")

# Set thresholds for significance
padj_threshold <- 0.05      # Adjusted p-value threshold (FDR)  
log2fc_threshold <- 1     # Log2 fold change threshold

cat(sprintf("Significance thresholds: FDR < %g, |logFC| > %g\n", 
            padj_threshold, log2fc_threshold))

# Get all significant genes (regardless of fold change direction)
sig_genes <- subset(results_df, FDR < padj_threshold & !is.na(FDR))

# Get upregulated genes (Treatment > Control)
up_genes <- subset(results_df, FDR < padj_threshold & 
                     logFC > log2fc_threshold & 
                     !is.na(FDR))

# Get downregulated genes (Treatment < Control)  
down_genes <- subset(results_df, FDR < padj_threshold & 
                       logFC < -log2fc_threshold & 
                       !is.na(FDR))

# Print summary statistics
cat(sprintf("Summary of results:\n"))
cat(sprintf("- Total significant genes: %d\n", nrow(sig_genes)))
cat(sprintf("- Upregulated genes: %d\n", nrow(up_genes)))
cat(sprintf("- Downregulated genes: %d\n", nrow(down_genes)))

# ==============================================================================
# 12. SAVE NORMALIZED EXPRESSION DATA
# ==============================================================================
cat("Saving normalized expression data...\n")

# Get log2 CPM values for normalized expression
normalized_expr <- cpm(dge_filtered, log = TRUE)

# Create dataframe with gene IDs
normalized_df <- as.data.frame(normalized_expr)
normalized_df$gene_id <- dge_filtered$genes$gene_id

# Reorder columns to put gene_id first
normalized_df <- normalized_df %>%
  select(gene_id, everything())

# Save normalized expression data
write.csv(normalized_df, file = "edgeR_normalized_logCPM.csv", row.names = FALSE)
cat("Normalized log2 CPM data saved to 'edgeR_normalized_logCPM.csv' ✓\n")

# Display first few rows for verification
cat("First few rows of normalized data:\n")
print(head(normalized_df, 3))

# ==============================================================================
# 13. FORMAT AND SAVE RESULTS
# ==============================================================================
cat("Formatting and saving results...\n")

# Function to format output dataframes consistently
format_output_df <- function(df, description) {
  cat(sprintf("Formatting %s (%d genes)...\n", description, nrow(df)))
  
  formatted_df <- df %>%
    select(gene_id, 
           logFC, 
           logCPM, 
           F, 
           PValue, 
           FDR) %>%
    arrange(FDR) %>%  # Sort by FDR (most significant first)
    mutate(across(where(is.numeric), ~round(.x, 6)))  # Round numeric values
  
  return(formatted_df)
}

# Format results
sig_genes_formatted <- format_output_df(sig_genes, "significant genes")
up_genes_formatted <- format_output_df(up_genes, "upregulated genes")  
down_genes_formatted <- format_output_df(down_genes, "downregulated genes")

# Save results to CSV files
write.csv(sig_genes_formatted, "edgeR_significant_genes.csv", row.names = FALSE)
write.csv(up_genes_formatted, "edgeR_upregulated_genes.csv", row.names = FALSE)
write.csv(down_genes_formatted, "edgeR_downregulated_genes.csv", row.names = FALSE)

# Save complete results table
write.csv(results_df, "edgeR_complete_results.csv", row.names = FALSE)

cat("Results saved:\n")
cat("- edgeR_significant_genes.csv\n")
cat("- edgeR_upregulated_genes.csv\n") 
cat("- edgeR_downregulated_genes.csv\n")
cat("- edgeR_complete_results.csv\n")

# ==============================================================================
# 14. CREATE VISUALIZATION PLOTS
# ==============================================================================
cat("Creating visualization plots...\n")

# 1. Volcano Plot
png("edgeR_Volcano_Plot_Treatment_vs_Control.png", width = 1000, height = 800, res = 300)

# Prepare data for volcano plot (convert to DESeq2-like format for EnhancedVolcano)
volcano_data <- results_df
volcano_data$padj <- volcano_data$FDR  # EnhancedVolcano expects 'padj'
volcano_data$log2FoldChange <- volcano_data$logFC  # EnhancedVolcano expects 'log2FoldChange'

volcano_plot <- EnhancedVolcano(volcano_data,
                                lab = volcano_data$gene_id,
                                x = 'log2FoldChange',
                                y = 'padj',
                                title = 'edgeR: Differential Gene Expression - Treatment vs Control',
                                subtitle = sprintf('Thresholds: FDR < %g, |logFC| > %g', 
                                                   padj_threshold, log2fc_threshold),
                                caption = sprintf('Total genes: %d | Significant: %d | Up: %d | Down: %d',
                                                  nrow(results_df), nrow(sig_genes), 
                                                  nrow(up_genes), nrow(down_genes)),
                                pCutoff = padj_threshold,
                                FCcutoff = log2fc_threshold,
                                pointSize = 2.0,
                                labSize = 3.0,
                                # Limit number of labels to avoid overcrowding
                                selectLab = c(head(up_genes_formatted$gene_id, 5), 
                                              head(down_genes_formatted$gene_id, 5)),
                                drawConnectors = TRUE,
                                widthConnectors = 0.5,
                                colConnectors = 'grey50')

print(volcano_plot)
dev.off()
cat("Volcano plot saved as 'edgeR_Volcano_Plot_Treatment_vs_Control.png' ✓\n")

# 2. MA Plot
png("edgeR_MA_plot.png", width = 800, height = 600, res = 300)
plotMD(qlf, main = "MA Plot: Treatment vs Control")
abline(h = c(-log2fc_threshold, log2fc_threshold), col = "red", lty = 2)
dev.off()
cat("MA plot saved as 'edgeR_MA_plot.png' ✓\n")

# 3. Heatmap of top differentially expressed genes
if (nrow(sig_genes) > 0) {
  cat("Creating heatmap of top significant genes...\n")
  
  # Get top 50 most significant genes (or all if fewer than 50)
  top_genes <- head(sig_genes_formatted$gene_id, min(50, nrow(sig_genes_formatted)))
  
  # Extract normalized expression for these genes
  heatmap_data <- normalized_expr[dge_filtered$genes$gene_id %in% top_genes, ]
  rownames(heatmap_data) <- dge_filtered$genes$gene_id[dge_filtered$genes$gene_id %in% top_genes]
  
  # Create annotation for samples
  annotation_col <- data.frame(
    Condition = col_data$condition,
    row.names = colnames(heatmap_data)
  )
  
  # Create heatmap
  png("edgeR_heatmap_top_genes.png", width = 1000, height = 800, res = 300)
  pheatmap(heatmap_data,
           annotation_col = annotation_col,
           scale = "row",
           clustering_distance_rows = "correlation",
           clustering_distance_cols = "correlation", 
           main = sprintf("Top %d Significant Genes", length(top_genes)),
           fontsize = 8,
           show_rownames = TRUE,
           show_colnames = TRUE)
  dev.off()
  cat("Heatmap saved as 'edgeR_heatmap_top_genes.png' ✓\n")
}

# ==============================================================================
# 15. ANALYSIS SUMMARY
# ==============================================================================
cat("\n" + paste(rep("=", 80), collapse = "") + "\n")
cat("edgeR ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat(paste(rep("=", 80), collapse = "") + "\n")
cat("Summary of outputs:\n")
cat(sprintf("- Normalized log2 CPM data: edgeR_normalized_logCPM.csv\n"))
cat(sprintf("- All significant genes (%d): edgeR_significant_genes.csv\n", nrow(sig_genes_formatted)))
cat(sprintf("- Upregulated genes (%d): edgeR_upregulated_genes.csv\n", nrow(up_genes_formatted)))
cat(sprintf("- Downregulated genes (%d): edgeR_downregulated_genes.csv\n", nrow(down_genes_formatted)))
cat(sprintf("- Complete results: edgeR_complete_results.csv\n"))
cat(sprintf("- Volcano plot: edgeR_Volcano_Plot_Treatment_vs_Control.png\n"))
cat(sprintf("- MA plot: edgeR_MA_plot.png\n"))
cat(sprintf("- MDS plot: MDS_plot.png\n"))
cat(sprintf("- BCV plot: BCV_plot.png\n"))
cat(sprintf("- Boxplot: Boxplot_logCPM.png\n"))
if (nrow(sig_genes) > 0) {
  cat(sprintf("- Heatmap: edgeR_heatmap_top_genes.png\n"))
}
cat(paste(rep("=", 80), collapse = "") + "\n")

# Display top 5 most significant genes
if (nrow(sig_genes_formatted) > 0) {
  cat("Top 5 most significant genes:\n")
  print(head(sig_genes_formatted, 5))
} else {
  cat("No significant genes found with current thresholds.\n")
  cat("Consider relaxing thresholds or checking data quality.\n")
}

# Display method comparison note
cat("\n" + paste(rep("-", 50), collapse = "") + "\n")
cat("NOTE: edgeR vs DESeq2 Differences:\n")
cat("- edgeR uses FDR instead of padj (same concept)\n")
cat("- edgeR uses logFC instead of log2FoldChange\n") 
cat("- edgeR uses TMM normalization vs DESeq2's median-of-ratios\n")
cat("- Results may differ slightly but should be broadly consistent\n")
cat(paste(rep("-", 50), collapse = "") + "\n")