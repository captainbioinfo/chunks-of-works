# Load required packages
library(ggplot2)
library(dplyr)
library(ggrepel)
# Load required packages
library(ggplot2)
library(dplyr)
library(ggrepel)

# Read in the data
data <- read.csv("top_100_pvalue.csv")
sampleinfo <- read.csv("sampleinfo.csv")

# Ensure columns are numeric (especially if read as characters)
data$log2FoldChange <- as.numeric(data$log2FoldChange)
data$pvalue <- as.numeric(as.character(data$pvalue))  # Convert to numeric, handling factors

# Replace NA values in pvalue with a large number (or any other strategy)
data$pvalue[is.na(data$pvalue)] <- 1

# Create the negLog10Pval column
data$negLog10Pval <- -log10(data$pvalue)

# Define thresholds for significance
log2fc_threshold <- 1
pvalue_threshold <- 0.05

# Create column for significance with categories
data$Significance <- with(data, case_when(
  abs(log2FoldChange) > log2fc_threshold & pvalue < pvalue_threshold & log2FoldChange > 0 ~ "Upregulated",
  abs(log2FoldChange) > log2fc_threshold & pvalue < pvalue_threshold & log2FoldChange < 0 ~ "Downregulated",
  TRUE ~ "Not Significant"
))

# Base plot
volcano <- ggplot(data, aes(x = log2FoldChange, y = negLog10Pval)) +
  geom_point(aes(color = Significance), alpha=0.7, size=2) +
  scale_color_manual(values=c("Upregulated" = "#e31a1c", "Downregulated" = "#1f78b4", "Not Significant" = "grey70")) +
  theme_minimal(base_size = 16, base_family = "Arial") +
  xlab(expression(Log[2]~Fold~Change)) +
  ylab(expression(-Log[10]~P~value)) +
  ggtitle("Volcano Plot: Differential Gene Expression") +
  theme(
    legend.title = element_blank(),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
    axis.title = element_text(face = "bold"),
    legend.position = "right"
  )

# Add labels for top significant genes
top_genes <- data[data$Significance != "Not Significant", ]
top_genes <- top_genes[order(top_genes$pvalue), ]
top_genes <- head(top_genes, 10)

volcano <- volcano +
  geom_text_repel(
    data=top_genes,
    aes(label=Gene_Name),
    size=4,
    box.padding = 0.5,
    point.padding = 0.5,
    segment.color = 'grey50',
    max.overlaps = Inf, # allow labeling more genes if overlaps occur
    seed = 42
  )

# Save plot as PNG at 600 dpi, size 8x6 inches (adjustable)
ggsave("volcano_plot_highres.png", volcano,
       width = 8, height = 6,
       dpi = 600, units = "in", bg = "white")

