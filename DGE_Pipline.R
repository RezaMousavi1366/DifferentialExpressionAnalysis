# Load required libraries
library(DESeq2)
library(tidyverse)

# Define helper function to reorder metadata
reorder_metadata <- function(metadata, count_data) {
  metadata <- metadata[match(colnames(count_data), rownames(metadata)), ]
  return(metadata)
}

# Define helper function for saving results
save_results <- function(results_data, output_path) {
  write.csv(as.data.frame(results_data), file = output_path, row.names = TRUE)
}

# Step 1: Data Preparation -----------------------------

# Load gene expression counts
gene_expression <- read.csv("D:/Reza/GenomicDataAnalysis/DESeq2DataAnlysis/counts_data.csv", row.names = 1)

# Load sample metadata
sample_metadata <- read.csv("D:/Reza/GenomicDataAnalysis/DESeq2DataAnlysis/sample_info.csv", row.names = 1)

# Validate column-row match
if (!all(colnames(gene_expression) %in% rownames(sample_metadata))) {
  stop("Column names in counts_data do not match row names in sample_info.")
}

# Reorder metadata to match gene expression columns
sample_metadata <- reorder_metadata(sample_metadata, gene_expression)

# Step 2: Construct DESeqDataSet -----------------------

dds <- DESeqDataSetFromMatrix(
  countData = gene_expression,
  colData = sample_metadata,
  design = ~ treatment
)

# Pre-filter low-count genes
dds <- dds[rowSums(counts(dds)) >= 10, ]

# Set untreated as the reference level for Treatment
dds$treatment <- relevel(dds$treatment, ref = "untreated")

# Step 3: Differential Expression Analysis -------------

# Run DESeq analysis
dds <- DESeq(dds)

# Extract results
results_default <- results(dds)

# Step 4: Explore and Save Results ---------------------

# Summarize the results
summary(results_default)

# Save the default results
save_results(results_default, "DESeq2_results_default.csv")

# Generate results with adjusted alpha
results_alpha <- results(dds, alpha = 0.01)
summary(results_alpha)
save_results(results_alpha, "DESeq2_results_alpha.csv")

# Perform contrasts
contrast_results <- results(dds, contrast = c("treatment", "treated_4hrs", "untreated"))
save_results(contrast_results, "DESeq2_contrast_results.csv")

# Step 5: Visualization --------------------------------

# Create MA plot
plotMA(results_default, ylim = c(-5, 5), main = "MA Plot of Differential Expression")

# Optional: Volcano Plot
volcano_data <- as.data.frame(results_default)
ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(aes(color = padj < 0.05), alpha = 0.8) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 Adjusted P-value")
