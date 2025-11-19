# Load the ggplot2 library
library(ggplot2)

# Create data frame
data <- data.frame(
  Dataset = rep(c("PRJNA544129", "PRJNA62742", "PRJNA43443"), each = 2),
  Category = rep(c("QC Passed (%)", "Mapped (%)"), times = 3),
  Percentage = c(81.74, 82.82, 99.77, 98.04, 99.44, 97.70)
)

# Convert Category to factor to control order and make fill mapping work
data$Category <- factor(data$Category, levels = c("QC Passed (%)", "Mapped (%)"))

# Plot
ggplot(data, aes(x = Dataset, y = Percentage, fill = Category)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f%%", Percentage)),
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3) +
  scale_fill_manual(values = c("QC Passed (%)" = "#E53935",  # red
                               "Mapped (%)" = "#1E88E5")) +  # blue
  labs(
    title = "Quality Control and Mapping Rates by Dataset",
    y = "Percentage (%)",
    x = "Dataset"
  ) +
  ylim(0, 110) +
  theme_minimal()
