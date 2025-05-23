# Load libraries
library(ggplot2)
library(DESeq2)
library(PCAtools)
library(readr)
library(ggrepel)

# Step 1: Read in gene expression data and sample metadata
counts_data <- read.csv("gene_count_matrix.csv", header = TRUE, sep = ",")  # Gene expression
col_data <- read.csv("sampleinfo.csv", header = TRUE, sep = ",")            # Sample metadata

# Step 2: Prepare count matrix
gene_ids <- make.unique(as.character(counts_data$gene_id))  # Ensure uniqueness
counts_data_no_gene_id <- counts_data[, -1]
rownames(counts_data_no_gene_id) <- gene_ids

# Step 3: Set rownames of metadata (must match colnames of count matrix)
rownames(col_data) <- col_data$gene_id  # ← adjust to your metadata column name
# Ensure column names in count data and rownames in metadata match
col_data <- col_data[colnames(counts_data_no_gene_id), ]

# Step 4: Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts_data_no_gene_id),
                              colData = col_data,
                              design = ~ condition)

# Step 5: Normalize using VST
dds <- DESeq(dds)
dds_vst <- vst(dds, blind = TRUE)
vsd_counts <- assay(dds_vst)

# Step 6: PCA on samples (transpose matrix)
pca_data <- prcomp(t(vsd_counts))  # ← sample-wise PCA

# Step 7: Variance explained
variance_explained <- round(100 * (pca_data$sdev^2 / sum(pca_data$sdev^2)), 2)

# Step 8: Prepare data for plotting
pca_plot_data <- data.frame(
  PC1 = pca_data$x[, 1],
  PC2 = pca_data$x[, 2],
  sample = rownames(pca_data$x),
  condition = col_data$condition[match(rownames(pca_data$x), rownames(col_data))]  # align metadata
)

# Create the PCA plot with merged color and shape legends
pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 5, alpha = 0.7) +  # Add transparency to points for better visibility
  geom_text_repel(aes(label = sample),  # Use ggrepel to avoid label overlap
                  max.overlaps = 100,   # Allow for more overlap of labels
                  size = 3,             # Adjust label size
                  box.padding = 0.5) +   # Add padding around the labels
  scale_color_manual(values = c("Control" = "blue", "Treatment" = "red")) +  # Color customization
  scale_shape_manual(values = c("Control" = 16, "Treatment" = 17)) +  # Control = circle, Treatment = triangle
  xlab(paste0("PC1 (", variance_explained[1], "% variance)")) +  # Label x-axis with variance explanation
  ylab(paste0("PC2 (", variance_explained[2], "% variance)")) +  # Label y-axis with variance explanation
  theme_minimal() +
  labs(title = "PRJNA627642", color = "Condition", shape = "Condition") +  # Combine color and shape legends into one
  theme(
    plot.title = element_text(hjust = 0.5,size = 24,  face = "bold"),  # Center the title
    legend.position = "top",  # Position legend to the top
    axis.text = element_text(size = 16),  # Adjust axis text size
    axis.title = element_text(size = 18),  # Adjust axis title size
    plot.background = element_rect(fill = "white", color = NA),  # Set background to white
    panel.background = element_rect(fill = "white", color = NA)  # Set panel background to white
  )

# Display the plot
print(pca_plot)
# Save the plot as a PNG for PowerPoint
#ggsave("pca_plot_for_ppt.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 600, bg = "white")

# Save the plot as a PDF for printing or high-quality documents
#ggsave("pca_plot_for_pdf22.png", plot = pca_plot, width = 3400, height = 3400, units = "px", dpi = 600, bg = "white")
################
# Add a new column for segment color based on condition
pca_plot_data$segment_color <- ifelse(pca_plot_data$condition == "Control", "blue", "red")

# Now create the plot
pca_plot <- ggplot(pca_plot_data, aes(x = PC1, y = PC2, color = condition, shape = condition)) +
  geom_point(size = 7, alpha = 0.7) +
  geom_text_repel(
    aes(label = sample, segment.color = segment_color),  # Use actual color names
    max.overlaps = 100,
    size = 5,
    box.padding = 0.5,
    segment.size = 0.5
  ) +
  scale_color_manual(values = c("Control" = "blue", "Treatment" = "red")) +
  scale_shape_manual(values = c("Control" = 16, "Treatment" = 17)) +
  xlab(paste0("PC1 (", variance_explained[1], "% variance)")) +
  ylab(paste0("PC2 (", variance_explained[2], "% variance)")) +
  theme_minimal() +
  labs(title = "PRJNA627642", color = "Condition", shape = "Condition") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 24, face = "bold"),
    legend.position = "top",
    legend.text = element_text(size = 20),
    legend.title = element_text(size = 22),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 22),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )
plot(pca_plot)
# Save the plot as a PNG for PowerPoint
ggsave("pca_plot_for_ppt.png", plot = pca_plot, width = 10, height = 8, units = "in", dpi = 600, bg = "white")

# Save the plot as a PDF for printing or high-quality documents
ggsave("pca_plot_for_pdf22.png", plot = pca_plot, width = 5400, height = 5400, units = "px", dpi = 600, bg = "white")
################
