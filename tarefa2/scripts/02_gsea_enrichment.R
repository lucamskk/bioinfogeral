library(DESeq2)
library(EnrichmentBrowser)
library(tidyverse)

# 1. Carregar Matriz de Contagens e Metadados
counts <- read.table("counts/count_matrix.tsv", header = TRUE, row.names = 1)
colData <- read.csv("colData.csv", row.names = 1)

colData$condition <- factor(colData$condition, levels = c("Control", "Treated"))
counts <- counts[, rownames(colData)]

# 2. Criar Objeto SummarizedExperiment
dds <- DESeqDataSetFromMatrix(countData = counts, colData = colData, design = ~ condition)
dds <- dds[rowSums(counts(dds)) >= 10, ]

se_data <- dds
colData(se_data)$GROUP <- ifelse(colData(se_data)$condition == "Treated", 1, 0)

# 3. Mapear Ensembl IDs para ENTREZID
se_data <- idMap(se_data, org = "hsa", from = "ENSEMBL", to = "ENTREZID")

# 4. Análise de Expressão Diferencial no SummarizedExperiment
se_data <- deAna(se_data, de.method = "DESeq2")

# 5. Baixar conjuntos de vias do KEGG para Homo sapiens
kegg_gs <- getGenesets(org = "hsa", db = "kegg")

# 6. Análise de Enriquecimento por Conjuntos (sbea com método gsea)
gsea_res <- sbea(method = "gsea", se = se_data, gs = kegg_gs, alpha = 0.05)

# 7. Salvar e Visualizar Resultados
gsea_table <- as.data.frame(gsea_res$res.tbl)
write.csv(gsea_table, "results/gsea_kegg_results.csv", row.names = FALSE)

# Exibir as 10 principais vias biológicas
head(gsea_table, 10)
