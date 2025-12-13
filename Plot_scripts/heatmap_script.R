# Load required libraries
library(ggplot2)
library(reshape2)
library(pheatmap)
library(RColorBrewer)

# Load data
data <- read.csv("all_c_gene_for_HP.csv", row.names = 1)

# Convert to numeric matrix
data_matrix <- as.matrix(data)
mode(data_matrix) <- "numeric"

# Log2 transformation to stabilize variance
data_log <- log2(data_matrix + 1)

# Identify Control vs Treatment from column names
group <- ifelse(grepl("_C", colnames(data)), "Control", "Treatment")
ann_col <- data.frame(Group = factor(group))
rownames(ann_col) <- colnames(data)

# Reorder columns to group Control and Treatment together
ordered_cols <- order(group)
data_log <- data_log[, ordered_cols]
ann_col <- ann_col[ordered_cols, , drop = FALSE]
group <- group[ordered_cols]

# Optionally, update column names with group label (cleaner visualization)
colnames(data_log) <- paste0(colnames(data_log), " (", group, ")")
rownames(ann_col) <- colnames(data_log)

# Use a softer, publication-quality color palette (e.g., RdBu from RColorBrewer)
heat_colors <- colorRampPalette(rev(brewer.pal(n = 11, name = "RdBu")))(100)
heat_colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(256)  # Increase resolution

# Annotation colors
ann_colors <- list(Group = c(Control = "#228B22", Treatment = "#B22222"))


# Plot and save heatmap as high-resolution PNG
png("heatmap22_shared_genes_cleaned.png", width = 7000, height = 3400, res = 600)
pheatmap(data_log,
         scale = "row",
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         color = heat_colors,
         annotation_col = ann_col,
         annotation_colors = ann_colors,
         fontsize_row = 12,
         fontsize_col = 8,
         fontsize = 16,
         border_color = NA,
         show_colnames = FALSE,
         main = "Expression Profile Heatmap of Common Genes")
dev.off()
########
