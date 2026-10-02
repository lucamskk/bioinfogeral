library(DESeq2)
library(tidyverse)
library(pheatmap)
library(ggrepel)
library(org.Hs.eg.db)

counts <- read.table("counts/count_matrix.tsv", header = TRUE, row.names = 1)
colData <- read.csv("colData.csv", row.names = 1)

colData$condition <- factor(colData$condition, levels = c("Control", "Treated"))
counts <- counts[, rownames(colData)]

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]
dds <- DESeq(dds)

res <- results(dds, contrast = c("condition", "Treated", "Control"))

df_res <- as.data.frame(res) %>%
    rownames_to_column("ensembl_id") %>%
    filter(!is.na(pvalue) & !is.na(log2FoldChange))

# Transforma em geneSymbol
df_res$symbol <- mapIds(org.Hs.eg.db, keys = df_res$ensembl_id, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
df_res <- df_res %>% mutate(gene_label = ifelse(is.na(symbol), ensembl_id, symbol))

write.csv(df_res, "results/deseq2_results.csv", row.names = FALSE)

# PCA
vsd <- vst(dds, blind = FALSE)
pca_plot <- plotPCA(vsd, intgroup = "condition") +
    theme_bw(base_size = 14) +
    labs(title = "PCA Plot - GSE296098 (Control vs Treated)") +
    theme(aspect.ratio = 1)
ggsave("results/pca_plot.png", pca_plot, width = 7, height = 7, dpi = 300)

# VolcanoPlot
df_res <- df_res %>% mutate(Status = case_when(
    log2FoldChange > 1 & padj < 0.05 ~ "Up-regulated",
    log2FoldChange < -1 & padj < 0.05 ~ "Down-regulated",
    TRUE ~ "Not Significant"
))

top_genes <- df_res %>% filter(Status != "Not Significant") %>% arrange(padj) %>% head(20)

volcano_gg <- ggplot(df_res, aes(x = log2FoldChange, y = -log10(padj), color = Status)) +
    geom_point(alpha = 0.5, size = 1.0) +
    scale_color_manual(values = c("Up-regulated" = "#d95f02", "Down-regulated" = "#7570b3", "Not Significant" = "grey80")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", alpha = 0.5) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
    geom_text_repel(data = top_genes, aes(label = gene_label), size = 3.5, fontface = "bold", max.overlaps = Inf) +
    theme_bw(base_size = 14) +
    labs(title = "Volcano Plot: Treated vs Control", x = expression(log[2] ~ "Fold Change"), y = expression(-log[10] ~ "(padj)"))
ggsave("results/volcano_plot.png", volcano_gg, width = 8, height = 8, dpi = 300)

# Heatmap top50 por padj
top50_degs <- df_res %>% filter(Status != "Not Significant") %>% arrange(padj) %>% head(50)
mat_vsd <- assay(vsd)[top50_degs$ensembl_id, ]
rownames(mat_vsd) <- top50_degs$gene_label

# normaliza por Z-score
pheatmap(mat_vsd, 
         scale = "row", 
         annotation_col = colData["condition"],
         show_colnames = TRUE,
         show_rownames = TRUE,
         main = "Heatmap - Top 50 DEGs (Z-score)",
         filename = "results/heatmap_degs.png",
         width = 8, height = 10)

