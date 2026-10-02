# Tarefa 2: Análise de Expressão Diferencial e GSEA em RNA-Seq (GSE296098)

Pipeline de bioinformática para processamento e análise estatística de RNA-Seq (*Homo sapiens*).

## Estrutura do Repositório
- `colData.csv`: Metadados das amostras (Control vs Treated - SRR33396499 a SRR33396504).
- `counts/count_matrix.tsv`: Matriz de contagens por gene gerada pelo `featureCounts`.
- `scripts/`: Scripts R para `DESeq2`, gráficos de controle e análise de enriquecimento de vias.
- `results/`: Gráficos gerados (`pca_plot.png`, `volcano_plot.png`).

# Apresentação por escrito já que não vai ser presencial ;(

## Parte 1: Escolha do dataset:
- Não escolhemos esse por nenhum motivo em específico, apenas porque sastisfez as requerições suas (bulk-rnaseq, 6 amostras)
- Através do FastQC/MultiQC verificamos que a qualidade por base estava excelente, com phred score > 30 (99,9% de precisão)
- A contaminação por adaptadores foi mínima (<0,1%) indicando que talvez as amostras já tenham passado por um controle de qualidade antes de terem sido depositadas no Banco de Dados
- Níveis de duplicação estavam na faixa de 30-45%, normais em análises de Bulk-RNASeq
- Alta profundidade de sequenciamento (> 50 Milhões de reads por amostra)
- %GC Content: 49-51%
- Nenhuma amostra se destoou das demais

## Parte 2: Alinhamento e contagem:
**Taxa de alinhamento** (HISAT2)
- **SRR33396499:** 96.97%
- **SRR33396500:** 97.81%
- **SRR33396501:** 97.55%
- **SRR33396502:** 95.89%
- **SRR33396503:** 95.91%
- **SRR33396504:** 95.39%

**featureCounts**
- **SRR33396499 (Control):** 68,6% (56.970.616 reads atribuídos)
- **SRR33396500 (Control):** 69,2% (55.324.180 reads atribuídos)
- **SRR33396501 (Control):** 69,1% (54.607.088 reads atribuídos)
- **SRR33396502 (Treated):** 65,9% (39.043.640 reads atribuídos)
- **SRR33396503 (Treated):** 65,9% (39.453.413 reads atribuídos)
- **SRR33396504 (Treated):** 65,5% (38.855.421 reads atribuídos)

## Parte 3: Expressão diferencial e GSEA no R
Utilizamos para expressão diferencial o DESeq2 e para GSEA o clusterProfiler. Esse útlimo é simplesmente por escolha pessoal minha (do lucas) pois já utilizo ele no dia a dia do laboratório :)
 - Através de uma leitura rápida do artigo original, nossos resultados parecem bater com os resultados achados pelos autores, encontrando em comum processos biológicos como regulação do metabolismo, reorganização da matriz extracelular e sinalização celular.
 - Vale salientar que através da análise GSEA (ranqueada) conseguimos encontrar nuances que não seriam encontradas em uma simples análise de sobreposição (ORA). Confesso que fiquei com vontade de testar uma análise topológica (usando netGSA) ou uma análise de rede (usando enrichmentBrowser), mas não tenho essa habilidade ainda kkkk

## Considerações finais:
Peço desculpas pela demora professor e agradeço muito a compreensão. Espero que confie que nós realmente tínhamos material para apresentar no dia correto, apenas não estavam no nível em que eu esperava (não batia muito com o artigo original, antes era outra GSE). Então como o senhor nos disponibilizou mais tempo para entregar o trabalho, decidimos recomeçar do zero com outros dados, dessa vez, dando um resultado muito mais satisfatório. Fiquei muito contente de finalmente realizar um trabalho utilizando técnicas que uso na Iniciação científica. Muito obrigado e de novo, peço desculpas pela demora, mas realmente a chuva complicou nossa vida.
