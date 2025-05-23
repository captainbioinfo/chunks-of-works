##################### Loading library #####
library(DESeq2)
library(tidyverse)
library(EnhancedVolcano)
######### Load count data and sample metadata #########################
counts_data <- read.csv("gene_count_matrix.csv")
col_data <- read.csv("pro_29_sampleinfo.csv")

# Store gene IDs before removing them from counts data
gene_ids <- counts_data$gene_id

# Remove gene_id column for DESeq2 processing
counts_data_no_gene_id <- counts_data[, -1]

# Ensure row names of count data match col_data
all(colnames(counts_data_no_gene_id) == rownames(col_data))

# Create DESeqDataSet
dds <- DESeqDataSetFromMatrix(countData = counts_data_no_gene_id,
                              colData = col_data,
                              design = ~ condition)

# Filter low-count genes
keep <- rowSums(counts(dds)) >= 73
dds <- dds[keep,]
# Also subset gene_ids to match filtered data
gene_ids <- gene_ids[keep]
# Relevel condition
dds$condition <- relevel(dds$condition, ref = "Control")
# Run DESeq2 analysis
dds <- DESeq(dds)
############ To save VST  ############################################################
# Perform VST transformation
vst_data <- vst(dds, blind = FALSE)
# Get the normalized expression data
vst_matrix <- assay(vst_data)
# Create a data frame with gene IDs and VST data
vst_df <- as.data.frame(vst_matrix)
vst_df$gene_id <- gene_ids
head_check <- head(vst_df, 3)
print("First few rows of VST data:")
print(head_check)
# Reorder columns to put gene_id first
vst_df <- vst_df %>%
  select(gene_id, everything())

# Save the VST-normalized data with gene IDs
write.csv(vst_df, file = "normalized_counts_data.csv", row.names = FALSE)

####################################################################################
res<- results(dds)
# Add gene IDs to results
res$gene_id <- gene_ids
#######################
##################### Loading libraries #####
library(DESeq2)
library(tidyverse)


############### Extract results for each condition vs. Control ###################
# Subset the DESeqDataSet to only include "Control" and "Treatment_24"
dds_subset <- dds[, dds$condition %in% c("Control", "Treatment")]
# Check 
table(dds_subset$condition)
####@@@@@@@@@@@@@@@@@@@

###############################################################################
#res <- results(dds)
res <- results(dds_subset, contrast = c("condition", "Control", "Treatment"))

# For a volcano plot:
res$gene_id <- gene_ids

#######
# Set thresholds
padj_threshold <- 0.05
log2fc_threshold <- 1.5
###
res_df <- as.data.frame(res)

res_df$gene_id <- gene_ids  # Use the gene_ids vector we stored earlier

################################################################################################
# Get significant genes
sig_genes <- subset(res_df, padj < padj_threshold)
up_genes <- subset(res_df, padj < padj_threshold & log2FoldChange > log2fc_threshold)
dim(up_genes); class(up_genes)


down_genes <- subset(res_df, padj < padj_threshold & log2FoldChange < -log2fc_threshold)
down_genes

####################
# Create function to format output data frame
format_output_df <- function(df) {
  df %>%
    select(gene_id, 
           baseMean, 
           log2FoldChange, 
           lfcSE, 
           stat, 
           pvalue, 
           padj) %>%
    arrange(padj)  # Sort by adjusted p-value
}

############################################
# Format all results
sig_genes_formatted <- format_output_df(sig_genes)
up_genes_formatted <- format_output_df(up_genes)
down_genes_formatted <- format_output_df(down_genes)

# Save results
write.csv(sig_genes_formatted, "significant.csv", row.names = FALSE)
write.csv(up_genes_formatted, "upregulated.csv", row.names = FALSE)
write.csv(down_genes_formatted, "downregulated.csv", row.names = FALSE)
##@@@@@@@@@@@@@@@@@

#########################################################################################################################
# Save volcano plot for Treatment_24 vs Control
png("Gene_Treatment_vs_Control.png")
EnhancedVolcano(res_df,
                lab = "", 
                x = 'log2FoldChange',
                y = 'padj',
                title = 'Treatment vs Control',
                pCutoff = padj_threshold,
                FCcutoff = log2fc_threshold,
                pointSize = 3.0, 
                labSize = 4.0)
dev.off()


##############################################################################################


