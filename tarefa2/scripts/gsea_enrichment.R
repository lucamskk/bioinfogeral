library(clusterProfiler)
library(enrichplot)
library(org.Hs.eg.db)
library(tidyverse)

df_res <- read.csv("results/deseq2_results.csv")

# ordena por log2FoldChange
df_ranked <- df_res %>%
    filter(!is.na(symbol) & !is.na(log2FoldChange)) %>%
    arrange(desc(log2FoldChange))

gene_list <- df_ranked$log2FoldChange
names(gene_list) <- df_ranked$symbol

# evita duplicado
gene_list <- gene_list[!duplicated(names(gene_list))]

# executa GSEA pelo GeneOntology
gsea_go <- gseGO(geneList     = gene_list,
                 OrgDb        = org.Hs.eg.db,
                 keyType      = "SYMBOL",
                 ont          = "BP", # determina qual grafo (Biological Process)
                 minGSSize    = 15,
                 maxGSSize    = 500,
                 pvalueCutoff = 0.05,
                 pAdjustMethod = "BH",
                 verbose      = FALSE)

# Dotplot
dotplot_fig <- dotplot(gsea_go, showCategory = 15, split = ".sign") + 
    facet_grid(.~.sign) +
    theme_bw(base_size = 12) +
    labs(title = "GSEA Dotplot: Vias Ativadas vs Suprimidas")
ggsave("results/gsea_dotplot.png", dotplot_fig, width = 10, height = 20, dpi = 300)

# gseaplot2 pelos top3 pathways
top3_pathways <- head(gsea_go@result$ID, 3)
gsea_top3 <- gseaplot2(gsea_go, geneSetID = top3_pathways, pvalue_table = TRUE)
ggsave("results/gsea_top3_pathways.png", gsea_top3, width = 11, height = 7, dpi = 300)
