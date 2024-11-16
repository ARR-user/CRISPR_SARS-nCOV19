# Load necessary libraries
library(DESeq2)
library(edgeR)
library(tidyverse)
library(pheatmap)
library(ggplot2)
library(ashr)


# Significance threshold for adjusted p-values (FDR)
signif_threshold <- 0.2

# Load data - counts and metadata
data <- read.delim("GSE159717_rnaseq_deseq_5dpi_counts_raw.tsv", header = TRUE, sep = "\t")

# Prepare count matrix (expression data)
count_data <- as.matrix(data[, grepl("^S_", colnames(data))])
rownames(count_data) <- data$gene_name

# Define conditions for the analysis (COVID-19 vs mock)
condition <- factor(c("Rem", "SARS", "mock", "Rem", "SARS", "mock"))
condition <- relevel(condition, ref = "SARS")  # Set "SARS" as the reference

# Create design matrix using the condition vector
design <- model.matrix(~ condition)
print(design)  # Check the design matrix


# ---------------------- edgeR Analysis ----------------------
# Create DGEList object for edgeR
dge <- DGEList(counts = count_data)

# Filter low-expression genes: keep genes with CPM > 1 in at least two samples
keep <- rowSums(cpm(dge) > 1) >= 2
dge <- dge[keep, , keep.lib.sizes = FALSE]  # Filtered DGEList object

# Normalize counts using TMM (to account for library size differences)
dge <- calcNormFactors(dge)

# Estimate dispersions
dge <- estimateDisp(dge, design)

# Fit generalized linear model (GLM) to normalized data
fit <- glmFit(dge, design)

# Likelihood ratio test for "mock" vs "SARS"
lrt_mock_vs_sars <- glmLRT(fit, coef = 2)

# Retrieve table of top genes by p-value
res_edgeR_mock_vs_sars <- topTags(lrt_mock_vs_sars, n = nrow(dge))$table

# Adjust p-values (Benjamini-Hochberg method)
res_edgeR_mock_vs_sars$FDR <- p.adjust(res_edgeR_mock_vs_sars$PValue, method = "BH")

# Filter for downregulated genes in SARS
downregulated_edgeR_mock_vs_sars <- rownames(res_edgeR_mock_vs_sars)[which(res_edgeR_mock_vs_sars$FDR < signif_threshold & res_edgeR_mock_vs_sars$logFC < 0)]
print(paste("Downregulated genes in edgeR (mock vs SARS):", length(downregulated_edgeR_mock_vs_sars)))


# ---------------------- DESeq2 Analysis ----------------------
# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = count_data, colData = data.frame(condition), design = ~ condition)

# Run DESeq2 pipeline (normalization, dispersion estimation, and differential expression)
dds <- DESeq(dds)

# Get results from DESeq2 for "mock" vs "SARS"
res_DESeq2_mock_vs_sars <- results(dds, contrast = c("condition", "mock", "SARS"), alpha = 0.3)

# Perform log fold change shrinkage (recommended for small sample sizes)
res_DESeq2_mock_vs_sars <- lfcShrink(dds, contrast = c("condition", "mock", "SARS"), res = res_DESeq2_mock_vs_sars, type = "ashr")
plotMA(res_DESeq2_mock_vs_sars, ylim=c(-10,10))

# Filter for downregulated genes in SARS
downregulated_DESeq2_mock_vs_sars <- rownames(res_DESeq2_mock_vs_sars)[which(res_DESeq2_mock_vs_sars$padj < signif_threshold & res_DESeq2_mock_vs_sars$log2FoldChange < 0)]
print(paste("Downregulated genes in DESeq2 (mock vs SARS):", length(downregulated_DESeq2_mock_vs_sars)))


# ---------------------- Combine Results ----------------------
# Combine results: Intersection of downregulated genes from both tools for mock vs SARS
common_downregulated_mock_vs_sars <- intersect(downregulated_DESeq2_mock_vs_sars, downregulated_edgeR_mock_vs_sars)

write.table(common_downregulated_mock_vs_sars, "common_downregulated_mock_vs_sars.txt", quote = FALSE, row.names = FALSE, col.names = FALSE)

print(paste("Common downregulated genes (mock vs SARS):", length(common_downregulated_mock_vs_sars)))


# ---------------------- Merge Results and Visualization ----------------------
# Pull log2 fold change and padj values for common downregulated genes in DESeq2
common_DESeq2_mock_vs_sars <- res_DESeq2_mock_vs_sars[common_downregulated_mock_vs_sars, c("log2FoldChange", "padj")]
common_DESeq2_mock_vs_sars <- as.data.frame(common_DESeq2_mock_vs_sars)
common_DESeq2_mock_vs_sars$Gene <- rownames(common_DESeq2_mock_vs_sars)

# Pull logFC and FDR values for common downregulated genes in edgeR
common_edgeR_mock_vs_sars <- res_edgeR_mock_vs_sars[common_downregulated_mock_vs_sars, c("logFC", "FDR")]
common_edgeR_mock_vs_sars <- as.data.frame(common_edgeR_mock_vs_sars)
common_edgeR_mock_vs_sars$Gene <- rownames(common_edgeR_mock_vs_sars)

# Merge results from both tools by gene
common_genes_merged_mock_vs_sars <- merge(common_DESeq2_mock_vs_sars, common_edgeR_mock_vs_sars, by = "Gene", all = TRUE)

# Sort and print the top downregulated genes
common_genes_sorted_mock_vs_sars <- common_genes_merged_mock_vs_sars[order(common_genes_merged_mock_vs_sars$padj), ]
print(head(common_genes_sorted_mock_vs_sars))


# ---------------------- PCA Plot for Sample Clustering ----------------------
# Perform PCA using variance stabilizing transformation (VST)
rld <- rlog(dds, blind = TRUE)

# PCA plot grouped by condition
plotPCA(rld, intgroup = "condition")


# ---------------------- Heatmap for Hierarchical Clustering ----------------------
# Compute sample distances based on the rlog-transformed data
sample_dists <- dist(t(assay(rld)))
sample_dist_matrix <- as.matrix(sample_dists)

annotation_col <- data.frame(condition = condition)
rownames(annotation_col) <- colnames(count_data)  # Assign sample names

# Heatmap of sample distances with hierarchical clustering
pheatmap(sample_dist_matrix, 
         clustering_distance_rows = sample_dists,
         clustering_distance_cols = sample_dists,
         annotation_col = annotation_col)
