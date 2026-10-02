# Tarefa 2: Análise de Expressão Diferencial e GSEA em RNA-Seq (GSE296098)

Pipeline de bioinformática para processamento e análise estatística de RNA-Seq (*Homo sapiens*).

## Estrutura do Repositório
- `colData.csv`: Metadados das amostras (Control vs Treated - SRR33396499 a SRR33396504).
- `counts/count_matrix.tsv`: Matriz de contagens por gene gerada pelo `featureCounts`.
- `scripts/`: Scripts R para `DESeq2`, gráficos de controle e análise de enriquecimento de vias.
- `results/`: Gráficos gerados (`pca_plot.png`, `volcano_plot.png`).
