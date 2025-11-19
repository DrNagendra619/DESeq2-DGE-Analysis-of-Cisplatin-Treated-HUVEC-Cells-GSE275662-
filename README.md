# DESeq2-DGE-Analysis-of-Cisplatin-Treated-HUVEC-Cells-GSE275662-
DESeq2 DGE Analysis of Cisplatin-Treated HUVEC Cells (GSE275662)
# 🧬 DESeq2 DGE Pipeline: Cisplatin Treatment in HUVEC Cells (GSE275662)

This R script automates a comprehensive bioinformatics workflow for **Differential Gene Expression (DGE)** analysis of the **GSE275662** RNA-Seq dataset. The study compares gene expression in Human Umbilical Vein Endothelial Cells (HUVEC) treated with **Cisplatin** against **Control** cells.

The pipeline utilizes the gold-standard **`DESeq2`** package, implements **`apeglm`** for accurate Log-Fold Change (LFC) estimation, and generates a suite of essential results and Quality Control (QC) visualizations.

## 🚀 Key Features

* **Automated Data Retrieval:** Downloads the gene count matrix directly from GEO using the provided file URL.
* **Robust DGE:** Implements the full `DESeq2` workflow, including pre-filtering and size factor normalization.
* **LFC Shrinkage:** Applies **`apeglm`** for more accurate and stabilized LFC estimates, especially for low-count genes.
* **Comprehensive Data Integration:** Correctly aligns, orders, and annotates the count matrix with custom metadata and gene names.
* **Integrated Visualization:** Generates three essential plots: **Volcano Plot**, **PCA Plot**, and **Heatmap** of the top DE genes.
* **Reproducible Output:** Saves the full, annotated DGE results to a CSV file.

---

## 🔬 Analysis Overview

| Component | Method / Test | Purpose |
| :--- | :--- | :--- |
| **Dataset** | GSE275662 | Effect of cisplatin on gene expression of HUVEC cells. |
| **DGE Tool** | `DESeq2` | Statistical method optimized for RNA-Seq count data. |
| **LFC Estimation** | `apeglm` | Provides robust and biologically meaningful LFCs. |
| **Comparison** | Cisplatin-Treated vs. Control | Identifies the transcriptional response of endothelial cells to the chemotherapeutic drug. |
| **QC Transformation** | **Regularized Log (rlog)** | Prepares data for global QC plots (PCA, Heatmaps). |

---

## 🛠️ Prerequisites and Setup

### 📦 Packages

The script automatically installs and loads the necessary Bioconductor and CRAN packages:
* `DESeq2`
* `SummarizedExperiment`
* `pheatmap`
* `EnhancedVolcano`
* `apeglm`
* `downloader`
* `tidyverse`
* `ggplot2`

### ⚙️ Execution

1.  **Download** the `DESeq2 DGE Analysis of Cisplatin-Treated HUVEC Cells (GSE275662).R` file.
2.  **Optional:** The output path is set to `D:/DOWNLOADS` (Step 1). You can change this path if needed.
3.  **Execute** the script in your R environment:
    ```R
    source("DESeq2 DGE Analysis of Cisplatin-Treated HUVEC Cells (GSE275662).R")
    ```

---

## 📁 Output Files (3 Plots + 1 CSV)

All output files are saved to the specified `output_path` (default: `D:/DOWNLOADS`).

### Statistical Results

| Filename | Type | Description |
| :--- | :--- | :--- |
| `DESeq2_results_GSE275662_Cisplatin.csv` | CSV | Full table of DGE results, including `gene_id`, `gene_name`, **shrunken log2FoldChange**, and **adjusted p-value** (padj). |

### Visualization and QC Plots

| Filename | Analysis Stage | Description |
| :--- | :--- | :--- |
| `Volcano_Plot_GSE275662.png` | Results | **Volcano Plot** showing shrunken $\log_2 \text{Fold Change}$ vs. $\log_{10}(\text{padj})$, highlighting DEGs. |
| `PCA_Plot_GSE275662.png` | QC / Results | **Principal Component Analysis (PCA)** plot demonstrating sample clustering and separation based on treatment condition. |
| `Heatmap_Top20_GSE275662.png` | Results | **Heatmap** visualizing the expression profiles of the **Top 20 Most Significant DE Genes** across all samples. |
