data <- read.csv("top_100_pvalue.csv")
sampleinfo <- read.csv("sampleinfo.csv")
head(data)
head(sampleinfo)

# Load required packages
library(pheatmap)
library(dplyr)
library(tibble)

# Read in the data
data <- read.csv("top_100_logfc.csv")
sampleinfo <- read.csv("sampleinfo.csv", header = FALSE)
colnames(sampleinfo) <- c("Sample", "Condition")

# Extract only normalized read counts columns
norm_counts <- data %>%
  select(Gene_Name, contains("Normalized.Read.Count"))

# Set gene names as row names
norm_counts <- column_to_rownames(norm_counts, var = "Gene_Name")

# Ensure the column order in expression matrix matches the sampleinfo
colnames(norm_counts) <- gsub("_Normalized.Read.Count", "", colnames(norm_counts))
sampleinfo$Sample <- gsub("_Raw.Read.Count", "", sampleinfo$Sample)
norm_counts <- norm_counts[, sampleinfo$Sample]

# Convert sample condition to row annotation for heatmap
annotation_col <- data.frame(Condition = sampleinfo$Condition)
rownames(annotation_col) <- sampleinfo$Sample

# Generate the heatmap
pheatmap(log2(norm_counts + 1),
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = TRUE,
         show_colnames = TRUE,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         fontsize_row = 6)
#####################
# Load libraries
library(pheatmap)
library(dplyr)
library(tibble)
library(RColorBrewer)

# Load data
data <- read.csv("top_100_pvalue.csv")
sampleinfo <- read.csv("sampleinfo.csv", header = FALSE)
colnames(sampleinfo) <- c("Sample", "Condition")

# Extract normalized counts
norm_counts <- data %>%
  select(Gene_Name, contains("Normalized.Read.Count"))

# Set gene names as rownames
norm_counts <- column_to_rownames(norm_counts, var = "Gene_Name")

# Clean and match sample names
colnames(norm_counts) <- gsub("_Normalized.Read.Count", "", colnames(norm_counts))
sampleinfo$Sample <- gsub("_Raw.Read.Count", "", sampleinfo$Sample)
norm_counts <- norm_counts[, sampleinfo$Sample]  # Reorder columns to match sample info

# Log2 transform normalized counts
log_counts <- log2(norm_counts + 1)

# Build annotation for columns
annotation_col <- data.frame(Group = sampleinfo$Condition)
rownames(annotation_col) <- sampleinfo$Sample

# Define color palette
ann_colors <- list(
  Group = c(control = "purple", T2D = "darkgreen")  # Dark Pink hex color
)
# Create the heatmap
heatmap <- pheatmap(log_counts,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "complete",
         border_color = NA,
         main = "Differentially Expressed Genes")
############%%%%%%%%
png("Heatmap_DEGs22.png", width = 8, height = 10, units = "in", res = 600)

pheatmap(log_counts,
         scale = "row",
         annotation_col = annotation_col,
         annotation_colors = ann_colors,
         color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 6,
         fontsize_col = 10,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         clustering_method = "complete",
         border_color = NA,
         main = "Differentially Expressed Genes")

dev.off()

