library(WGCNA)
library(flashClust)
library(ggplot2)
library(DESeq2)
library(ggplot2)
library(gridExtra)  # Needed for grid.arrange()
library(tidyverse)
library(gridExtra)
library(magrittr)
rm(list=ls())
enableWGCNAThreads()

data <- readr::read_delim("for_wgcna_project.txt",  delim = "\t")    # <= path to the data file
# ########To remove gene_id2 column #######################
 
sampleinfo <- read.table("sampleinfo.txt", header=T,sep="\t")    # <= path to the data file
data[1:5,1:10] ;str(data)                          
names(data)[1] = "GeneId"
names(data)           # Look at the column names

# ======== R processing to clean and tidy the dataset for exploratory graphics.

col_sel = names(data)[-1]     # Get all but first column nam
mdata <- data %>%
  tidyr::pivot_longer(
    .,                        # The dot is the the input data, magrittr tutorial
    col = all_of(col_sel)
  ) %>%  
  mutate(
    group = gsub("-.*","", name) %>% gsub("[.].*","", .)   # Get the shorter treatment names
  )

#========== Normalisation

de_input = as.matrix(data[,-1])
row.names(de_input) = data$GeneId
de_input[1:5,1:10]

meta_df <- data.frame( Sample = names(data[-1])) %>%
  mutate(
    Type = gsub("-.*","", Sample) %>% gsub("[.].*","", .)
  )
meta_df<-sampleinfo
dds <- DESeqDataSetFromMatrix(round(de_input),
                              meta_df,
                              design = ~ condition)##i change it from TYPE 

#dds <- DESeq(dds)
dds <- DESeq(dds, fitType = "local")

vsd <- varianceStabilizingTransformation(dds)
wpn_vsd <- getVarianceStabilizedData(dds)
rv_wpn <- rowVars(wpn_vsd)
summary(rv_wpn)
q75_wpn <- quantile( rowVars(wpn_vsd), .75)  # <= original
#q95_wpn <- quantile( rowVars(wpn_vsd), .95)  # <= changed to 95 quantile to reduce dataset
expr_normalized <- wpn_vsd[ rv_wpn > q75_wpn, ];
dim(expr_normalized)

#==================== WGCNA =================
input_mat = t(expr_normalized)
# Choose a set of soft-thresholding powers
#############
sampleTree = hclust(dist(input_mat), method = "average");
plot(sampleTree)
##################
# @@@@@@@@@@@@@@@@@@@@@@@@@ FOR PLOTTING OF SAMPLE TREE DENDOGRAM #@@@@@@@@@@@@@@@@@@@
# Perform hierarchical clustering (if not already done)
# Close any open graphics devices
while (dev.cur() > 1) dev.off()

# Recalculate clustering if needed
sampleTree <- hclust(dist(input_mat), method = "average")

# Save as wider PNG for clarity
png("sample_clustering_tree_clean.png", width = 3400, height = 2400, res = 300)

# Plot with more readable parameters
plot(sampleTree,
     main = "Sample Clustering Dendrogram",
     xlab = "",
     sub = "",
     cex = 0.9)               # Shrink label size slightly)

dev.off()

############## @@@@@@@@@@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@
powers = c(c(1:10), seq(from = 2, to=20, by=2));

#powers = c(5:30)  # Expanding the range of powers
# Calculate a soft-thresholding power used for downstream analysis
sft = pickSoftThreshold(input_mat, powerVector = powers, verbose = 5)
sft[["powerEstimate"]]

# Plot Scale-free topology fit index as a function of the soft-thresholding power
ggplot(sft$fitIndices, aes(Power, SFT.R.sq, label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  geom_hline(yintercept = 0.8, color = 'red') +
  labs(x = 'Power', y = 'Scale free topology model fit, signed R^2') +
  theme_classic()

# Plot Mean connectivity as a function of the soft-thresholding power
ggplot(sft$fitIndices, aes(Power, mean.k., label = Power)) +
  geom_point() +
  geom_text(nudge_y = 0.1) +
  labs(x = 'Power', y = 'Mean Connectivity') +
  theme_classic()

#@@@@@@@@@@@@@@@@@@@@@@@@@@################

#=============== WGCNA : Network formation

picked_power2 = sft$powerEstimate
print(picked_power2)
picked_power <- 6
temp_cor <- cor      
cor <- WGCNA::cor         # Force it to use WGCNA cor function (fix a namespace conflict issue)
netwk <- blockwiseModules(input_mat,                # <= input here
                          
                          # == Adjacency Function ==
                          power = picked_power,                # <= power here
                          networkType = "unsigned",
                          
                          # == Tree and Block Options ==
                          deepSplit = 2,
                          pamRespectsDendro = F,
                          # detectCutHeight = 0.75,
                          minModuleSize = 20, ## sir changed it from 30 to 20
                          maxBlockSize = 4000,
                          
                          # == Module Adjustments ==
                          reassignThreshold = 0,
                          mergeCutHeight = 0.25,
                          
                          # == TOM == Archive the run results in TOM file (saves time)
                          saveTOMs = T,
                          saveTOMFileBase = "FPKM-TOM",
                          
                 
                          # == Output Options
                          numericLabels = T,
                          verbose = 3)
cor <- temp_cor
##
#
moduleColors = labels2colors(netwk$colors)
#
MEs0 = moduleEigengenes(input_mat, moduleColors)$eigengenes
MEs = orderMEs(MEs0); #

# 
#file.remove('All_Gene_KME.txt')
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=as.data.frame(MEs[,paste("ME",module,sep="")])
  colnames(ME)=module
  datModExpr=input_mat[,moduleColors==module]
  datKME = signedKME(datModExpr, ME)
  datKME=cbind(datKME,rep(module,length(datKME)))
  write.table(datKME,quote = F,row.names = T,append = T,file = "All_Gene_KME.txt",col.names = F)
}

#===============模块绘图==================
expColor=t(numbers2colors(log10(input_mat+1),colors=blueWhiteRed(100),naColor="grey"))
colnames(expColor)=rownames(input_mat)
################## @@@@@ BETTER QULAITY PLOTS @@@@ #################################
######### PLOTING OF CLUSTER DENDOGRAM WITH HEATMAP #####################
# Save high-res plot (e.g., for publication or presentations)
png("wgcna.dendroColors.highres.png", height = 9000, width = 7400, res = 600)

# Plot with improved aesthetics
plotDendroAndColors(
  netwk$dendrograms[[1]],
  colors = cbind(moduleColors[netwk$blockGenes[[1]]], expColor),
  groupLabels = c("Module", colnames(expColor)),
  dendroLabels = FALSE,
  hang = 0.03,
  addGuide = TRUE,
  guideHang = 0.05,
  cex.rowText = 1.0,        # Bigger row label text
  main = "Gene Dendrogram and Module Colors"
)
dev.off()
###@@@@ PLOTING OF EIGENGINE ADJACENCY HEATMAP @@@###############################
png("wgcna.eigengene_adjacency.heatmapres.png", width = 5600, height = 5600, res = 600)
# Plot eigengene adjacency heatmap
plotEigengeneNetworks(
  MEs,
  setLabels = "Eigengene Adjacency Heatmap",
  plotDendrograms = FALSE,
  marHeatmap = c(8, 8, 4, 4),       # More space for axis labels
  cex.lab = 1.5,                    # Increase label size
  cex.main = 1.5                    # Increase title size
)
dev.off()
################ @@@@@@@@@@@@@@@@@@@@ @@@@@@@@@@@@@@@@@@@@@@@@
nSamples = ncol(de_input)
###############
for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  ME=MEs[,paste("ME",module,sep="")]
  png(paste("wgcna.", module, ".express.barplot.png", sep=""),height = 700,width = 900)
  par(mfrow=c(2,1),mar=c(0.3,5.5,3,2))
  plotMat(t(scale(input_mat[,moduleColors==module])),
          rlabels=F,main=module,cex.main=2,clabels=F)
  
  par(mar=c(5,4.2,0,0.7))
  barplot(ME,col=module,main="",cex.main=2,ylab="eigengene expression",xlab="sample")
  dev.off()
}
load('FPKM-TOM-block.1.RData') 
TOM=as.matrix(TOM)
#
TOM = TOMsimilarityFromExpr(input_mat, power =picked_power,TOMType = "unsigned"); 

for(module in substring(colnames(MEs),3)){
  if(module == "grey") next
  probes = colnames(input_mat)
  inModule = is.finite(match(moduleColors, module))
  modProbes = probes[inModule]
  modTOM = TOM[inModule, inModule]
  dimnames(modTOM) = list(modProbes, modProbes)
  cyt = exportNetworkToCytoscape(modTOM,
                                 edgeFile = paste("New folder/CytoscapeInput5-edges-", module, ".txt", sep=""),
                                 nodeFile = paste("New folder/CytoscapeInput5-nodes-", module, ".txt", sep=""),
                                 weighted = TRUE,
                                 threshold = 0.25,
                                 nodeNames = modProbes,
                                 nodeAttr = moduleColors[inModule])
}
###########################
allTraits<- read.table("pro_42_sampleinfo.txt", header=T,sep="\t")    # <= path to the data file
# Get sample names from your dataset

# First, let's examine your data structure
head(input_mat) # Check sample names in your expression matrix
head(allTraits) # As you showed, this contains SampleName, condition, time

# The problems are:
# 1. Matching samples between input_mat and allTraits
# 2. Converting non-numeric trait data to usable format

# Get sample names from your dataset
SampleName <- rownames(input_mat)

# Check if SampleName matches with allTraits$SampleName
match_result <- match(SampleName, allTraits$SampleName)
print(match_result) # Check if we have NA values or duplicates

# Alternative approach - make condition and time as factors and convert to numeric
# First, create proper data frame
datTraits <- data.frame(
  sample = SampleName,
  stringsAsFactors = FALSE
)

# Match with allTraits
for(i in 1:nrow(datTraits)) {
  idx <- which(allTraits$SampleName == SampleName[i])
  if(length(idx) > 0) {
    datTraits$condition[i] <- allTraits$condition[idx[1]]
    datTraits$time[i] <- allTraits$time[idx[1]]
  }
}

# Convert categorical variables to numeric for correlation
datTraits$condition_num <- as.numeric(as.factor(datTraits$condition))
#datTraits$time_num <- as.numeric(as.factor(datTraits$time))

# Use only numeric columns for correlation
#datTraits_numeric <- datTraits[, c("condition_num", "time_num")]
datTraits_numeric <- datTraits["condition_num"]

rownames(datTraits_numeric) <- SampleName

# Now calculate correlations
moduleTraitCor <- cor(MEs, datTraits_numeric, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

# Format text for the heatmap cells
textMatrix <- paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

# Create heatmap with better column labels
trait_labels <- c("Condition") # Better labels for the heatmap


###################### @@@@@@@@@@@ 
# Save the heatmap with high resolution
png("wgcna.Module-traitres.heatmap.png", width = 4000, height = 6000, res = 600)
# Set margins and plot the heatmap
par(mar = c(6, 10, 4, 2))  # Increased bottom and left margins for label space
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = trait_labels,
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(100),       # Smoother gradient
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 1.2,                   # Increased text size
  zlim = c(-1, 1),
  main = "Module-Trait Relationships"
)

dev.off()
################## @@@@@@@@@@@@@
# Save the results
write.csv(moduleTraitCor, "ModuleTraitCorrelations.csv")
write.csv(moduleTraitPvalue, "ModuleTraitPvalues.csv")


################################
# Create trait data frame properly matched to input_mat samples
# Get sample names from your expression data
sample_names <- rownames(input_mat)

# Create a new trait data frame with matching samples
trait_data <- data.frame(row.names = sample_names)

# Add columns for condition and time by matching with allTraits
trait_data$condition_num <- NA
trait_data$time_num <- NA

# Map values from allTraits to trait_data
for(i in 1:nrow(trait_data)) {
  # Find this sample in allTraits
  idx <- which(allTraits$SampleName == rownames(trait_data)[i])
  if(length(idx) > 0) {
    # Map condition and time
    trait_data$condition_num[i] <- as.numeric(as.factor(allTraits$condition))[idx[1]]
    trait_data$time_num[i] <- as.numeric(as.factor(allTraits$time))[idx[1]]
  }
}

# Check if we have any missing values
cat("Missing condition values:", sum(is.na(trait_data$condition_num)), "\n")
cat("Missing time values:", sum(is.na(trait_data$time_num)), "\n")

# If there are missing values, we might need to handle them
# For now, let's impute with the mean (or you can remove those samples)
trait_data$condition_num[is.na(trait_data$condition_num)] <- mean(trait_data$condition_num, na.rm=TRUE)
trait_data$time_num[is.na(trait_data$time_num)] <- mean(trait_data$time_num, na.rm=TRUE)

# Now calculate Gene Significance
trait_of_interest <- "condition_num"  # or "time_num"

#################################
# Get the trait vectors
trait_vector_condition <- as.numeric(trait_data[, "condition_num"])
trait_vector_time <- as.numeric(trait_data[, "time_num"])

# Calculate gene significance for condition
gene_significance_condition <- numeric(ncol(input_mat))
for(i in 1:ncol(input_mat)) {
  gene_significance_condition[i] <- abs(cor(input_mat[, i], trait_vector_condition, use="p"))
}

# Calculate gene significance for time
gene_significance_time <- numeric(ncol(input_mat))
for(i in 1:ncol(input_mat)) {
  gene_significance_time[i] <- abs(cor(input_mat[, i], trait_vector_time, use="p"))
}

# Create data frames with gene module colors and gene significance
module_GS_condition <- data.frame(
  GS = gene_significance_condition,
  module = moduleColors
)

module_GS_time <- data.frame(
  GS = gene_significance_time,
  module = moduleColors
)

# Calculate the mean Gene Significance for each module
mean_GS_condition <- tapply(module_GS_condition$GS, module_GS_condition$module, mean)
mean_GS_time <- tapply(module_GS_time$GS, module_GS_time$module, mean)

# Get unique module colors (excluding grey)
modules <- names(mean_GS_condition)
if("grey" %in% modules) {
  modules <- modules[modules != "grey"]
}

# Create a data frame for plotting
plot_data <- data.frame(
  Module = modules,
  Condition = mean_GS_condition[modules],
  Time = mean_GS_time[modules]
)
########################## GENE SIGNIFICANCE MODULE IN CONDITION TRAIT ONLY ################################################################################################
##################
# Set up high-resolution PNG output
png("wgcna2.combined_gene_significance_condition.png", width = 5500, height = 4400, res = 600)

# Set up a single plotting layout
layout(matrix(1, nrow = 1))

# Adjust margins: bottom, left, top, right
par(mar = c(8, 7, 5, 2))

# Calculate max y-value for y-axis scaling
y_max <- max(plot_data$Condition, na.rm = TRUE) * 1.2

# Create the barplot
bp_pos <- barplot(plot_data$Condition,
                  col = "darkblue",
                  border = "white",
                  names.arg = rep("", length(plot_data$Module)),  # Suppress x labels here
                  ylim = c(0, y_max),
                  main = "Mean Gene Significance by Module (Condition)",
                  ylab = "Mean Gene Significance",
                  cex.main = 2.2,          # Main title size
                  cex.lab = 1.8,           # Axis label size
                  cex.axis = 1.5           # Axis tick size
)

# Add module names rotated below bars
text(x = bp_pos, y = par("usr")[3] - 0.02 * y_max,
     labels = plot_data$Module,
     srt = 45, adj = 1, xpd = TRUE,
     cex = 1.5, font = 2, col = "black")


# Finish
dev.off()
######@@@@@@@@@@@
# High-res output
png("wgcna2kjglg.combined_gene_significance_condition.png", width = 5500, height = 4400, res = 600)

layout(matrix(1, nrow = 1))
par(mar = c(8, 7, 5, 2))

y_max <- max(plot_data$Condition, na.rm = TRUE) * 1.2

# Barplot with WGCNA module colors
bp_pos <- barplot(plot_data$Condition,
                  col = plot_data$Module,    # Module color as bar fill
                  border = "white",
                  names.arg = rep("", length(plot_data$Module)),
                  ylim = c(0, y_max),
                  main = "Mean Gene Significance by Module (Condition)",
                  ylab = "Mean Gene Significance",
                  cex.main = 2.2,
                  cex.lab = 1.8,
                  cex.axis = 1.5)

# Add rotated module names below bars
text(x = bp_pos, y = par("usr")[3] - 0.02 * y_max,
     labels = plot_data$Module,
     srt = 45, adj = 1, xpd = TRUE,
     cex = 1.5, font = 2, col = "black")

dev.off()

