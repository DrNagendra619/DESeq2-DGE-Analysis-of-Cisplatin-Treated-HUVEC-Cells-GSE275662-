### -------------------------------------------
### FINAL SCRIPT: DESeq2 DGE Analysis (Cisplatin / GSE275662)
### -------------------------------------------
###
### TITLE: Effect of cisplatin on gene expression of human umbilical vein endothelial cells
### DATASET: GSE275662
###
### -------------------------------------------

## -------------------------------------------
## Step 0: Install and Load Packages
## -------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(
  c("DESeq2", "SummarizedExperiment", "pheatmap", "EnhancedVolcano", "apeglm", "downloader"),
  update = FALSE, ask = FALSE
)

# Install tidyverse for data manipulation
install.packages("tidyverse")

# Load libraries
library(DESeq2)
library(SummarizedExperiment)
library(tidyverse)        # For data cleaning and plotting
library(downloader)       # For downloading the data file
library(pheatmap)         # For heatmaps
library(EnhancedVolcano)  # For volcano plots
library(apeglm)           # For LFC shrinkage
library(ggplot2)          # For plotting

## -------------------------------------------
## Step 1: CONFIGURATION - Set Output Path
## -------------------------------------------
# Define the default folder where all results will be saved.
# We use forward slashes "/" for R paths, even on Windows.
output_path <- "D:/DOWNLOADS"

# Create the directory if it doesn't already exist
if (!dir.exists(output_path)) {
  dir.create(output_path, recursive = TRUE)
  print(paste("Created output directory:", output_path))
} else {
  print(paste("Output directory already exists:", output_path))
}

## -------------------------------------------
## Step 2: Download and Load Data (CORRECTED)
## -------------------------------------------
# URL for the gene count matrix from GEO
data_url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE275662&format=file&file=GSE275662%5Fgene%2Ecount%2Ematrix%2Eannot%2Etxt%2Egz"
data_file <- "GSE275662_gene.count.matrix.annot.txt.gz"

# Download the file if it doesn't exist
if (!file.exists(data_file)) {
  download(data_url, destfile = data_file, mode = "wb")
}

# Load the count data
# We now know the file is a simple t-separated file (tsv)
# The first column ('gene_id') will be used as row names
count_data_annotated <- read.delim(data_file,
                                   header = TRUE,
                                   row.names = 1,
                                   check.names = FALSE,
                                   sep = "\t") # Specify tab-separated

## -------------------------------------------
## Step 3: Prepare Count Matrix and Metadata (CORRECTED)
## -------------------------------------------

# FIRST: Create the sample metadata (colData)
# We use the *ACTUAL* names from the file header.
# Note the order is mixed, so we must match it.
sample_info <- data.frame(
  condition = factor(c("cisplatin", "control", "control", "control", "cisplatin", "cisplatin")),
  row.names = c(
    "HUVEC_1_3", # Cisplatin
    "HUVEC_0_2", # Control
    "HUVEC_0_1", # Control
    "HUVEC_0_3", # Control
    "HUVEC_1_1", # Cisplatin
    "HUVEC_1_2"  # Cisplatin
  )
)

# Get the list of sample names from our metadata
sample_names <- rownames(sample_info)

# SECOND: Create the count_matrix
# This will now work because `sample_names` matches the file header
count_matrix <- count_data_annotated %>%
  select(all_of(sample_names))

# Make sure the counts are integers
count_matrix <- round(count_matrix)

# THIRD: Verify and Re-order
# Ensure count_matrix columns are in the *exact same order* as sample_info rows
if (!all(colnames(count_matrix) == rownames(sample_info))) {
  # Re-order count_matrix columns to match sample_info rows
  count_matrix <- count_matrix[, rownames(sample_info)]
}

# Final verification
if (!all(colnames(count_matrix) == rownames(sample_info))) {
  stop("Error: Column names in count_matrix still do not match row names in sample_info.")
}

print("✅ Step 3 complete: count_matrix and sample_info are perfectly matched.")

## -------------------------------------------
## Step 4: Create DESeq2 Object
## -------------------------------------------
# This step will now work
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = sample_info,
  design = ~ condition
)

# Pre-filter rows with very few reads (removes genes with < 10 total reads)
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep,]

## -------------------------------------------
## Step 5: Run DESeq2 Analysis (CORRECTED)
## -------------------------------------------
# Run the main DESeq2 function
dds <- DESeq(dds)

# Check the names of the comparisons (coefficients)
# This will now be "condition_control_vs_cisplatin" because we set the ref
print(resultsNames(dds))

# ***THIS IS THE FIX***:
# We set coef_name to *exactly* match the output from resultsNames()
coef_name <- "condition_control_vs_cisplatin"

# Apply Log-Fold Change (LFC) shrinkage
# This will now work
resLFC <- lfcShrink(dds, coef = coef_name, type = "apeglm")


## -------------------------------------------
## Step 6: Annotate and Save Results
## -------------------------------------------
# Get the gene names from the original downloaded file
# We select the 'gene_name' column
gene_annotations <- count_data_annotated %>%
  select(gene_name) %>%
  rownames_to_column(var = "gene_id")

# Convert LFC results to a data frame
res_data <- as.data.frame(resLFC)
res_data <- rownames_to_column(res_data, var = "gene_id")

# Join the results with the gene name annotations
final_results <- left_join(res_data, gene_annotations, by = "gene_id")

# Order the final results by adjusted p-value
final_results_ordered <- final_results[order(final_results$padj), ]

# Show top 10 results
print("--- Top 10 Differentially Expressed Genes ---")
final_results_ordered %>%
  select(gene_id, gene_name, log2FoldChange, pvalue, padj) %>%
  head(10)

# Save the full results table to a CSV file
results_csv_file <- file.path(output_path, "DESeq2_results_GSE275662_Cisplatin.csv")
write.csv(final_results_ordered, results_csv_file, row.names = FALSE)

## -------------------------------------------
## Step 7: Volcano Plot (and Save)
## -------------------------------------------
# Create a volcano plot
volcano_plot <- EnhancedVolcano(
  final_results_ordered,
  lab = final_results_ordered$gene_name, # Use gene_name for labels
  x = "log2FoldChange",
  y = "padj", # Use adjusted p-value
  pCutoff = 0.05,
  FCcutoff = 1,
  title = "Volcano Plot: Cisplatin vs Control (GSE275662)"
)

# Save the plot
ggsave(
  filename = file.path(output_path, "Volcano_Plot_GSE275662.png"),
  plot = volcano_plot,
  width = 10,
  height = 8
)

## -------------------------------------------
## Step 8: PCA Plot (and Save)
## -------------------------------------------
# Transform the count data for visualization
rld <- rlog(dds)

# Generate a Principal Component Analysis (PCA) plot
pca_plot <- plotPCA(rld, intgroup = "condition") +
  theme_bw() +
  ggtitle("PCA Plot: GSE275662")

# Save the plot
ggsave(
  filename = file.path(output_path, "PCA_Plot_GSE275662.png"),
  plot = pca_plot,
  width = 8,
  height = 6
)

## -------------------------------------------
## Step 9: Heatmap of Top 20 DE Genes (and Save)
## -------------------------------------------
# Get the gene IDs of the top 20 most significant genes
top20_genes <- final_results_ordered$gene_id[1:20]

# Get the transformed counts for these genes
mat <- assay(rld)[top20_genes, ]

# Z-score scaling: Scale by row (gene)
mat <- mat - rowMeans(mat)

# Define the file path for the heatmap
heatmap_file <- file.path(output_path, "Heatmap_Top20_GSE275662.png")

# Create a heatmap and save it to the file
pheatmap(
  mat,
  scale = "row", # Already scaled, but good practice
  main = "Top 20 Differentially Expressed Genes (GSE275662)",
  filename = heatmap_file,
  width = 8,
  height = 10
)

## -------------------------------------------
## Step 10: Completion
## -------------------------------------------
message("---")
message(paste("✅ Done! All results and plots saved to:", output_path))
message("---")
