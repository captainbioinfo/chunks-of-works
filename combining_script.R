# ==============================================================================
# Merge multiple count files column-wise using GeneId
# ==============================================================================

# Load necessary packages
# install.packages("dplyr") # if not installed
library(dplyr)

# Step 1: Read all CSVs
f1 <- read.csv("PRJNA544129_edited.csv")
f2 <- read.csv("PRJEN43443_edited.csv")
f3 <- read.csv("PRJNA551141_edited.csv")
f4 <- read.csv("PRJNA556769_edited.csv")
f5 <- read.csv("PRJNA591729_edited.csv")
f6 <- read.csv("PRJNA627642_edited.csv")
f7 <- read.csv("PRJNA778892_edited.csv")

# Step 2: Merge all data frames by "GeneId"
# Use reduce + merge with 'by = "GeneId"' and all = TRUE to preserve all GeneIds
all_files <- list(f1, f2, f3, f4, f5, f6, f7)

merged_data <- Reduce(function(x, y) merge(x, y, by = "GeneId", all = TRUE), all_files)

# Step 3: Replace NA with 0
merged_data[is.na(merged_data)] <- 0

# Step 4: Preview
cat("✅ Merged table preview:\n")
print(head(merged_data))

# Step 5: Save output
write.csv(merged_data, "Merged_By_GeneId_AllSamples.csv", row.names = FALSE)
cat("✔ File saved: Merged_By_GeneId_AllSamples.csv\n")
