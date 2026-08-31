# Tarefa 1 para a cadeira de Bioinformática Geral

## Objetivo:
- **Organismo:** *Cupriavidus taiwanensis* LMG 19424
- **Estrutura genômica:** DNA
- **Assembly/release de referência:** GCA_000069785.1 / GCF_000069785.1, assembly `ASM6978v1` (NCBI); disponível também no Ensembl Bacteria como `Cupriavidus_taiwanensis_lmg_19424_gca_000069785`.
- **Justificativa da escolha:** Aleatório
- **Como os dados foram obtidos:**  https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/fasta/bacteria_105_collection/cupriavidus_taiwanensis_gca_900249975/dna/Cupriavidus_taiwanensis_gca_900249975.CBM2634_.dna.toplevel.fa.gz

## Requisitos e dependências
 
Ferramentas utilizadas neste pipeline (ajuste as versões conforme o que você usar):
 
| Ferramenta | Versão | Finalidade |
|---|---|---|
| `sra-tools`| `3.4.1` | Download dos reads |
| `bwa` | `0.7.19` | Alinhamento |
| `samtools` | `1.24` | Indexação/manipulação de BAM |

##Download das READS

> Fonte: SRR3943818
```bash
mkdir -p reads
cd reads
sudo apt update && sudo apt install -y sra-toolkit
fastq-dump -X 100000 --split-files SRR3943818
```
##Download do Genoma de Referência

> Fonte: Ensembl Bacteria
> Assembly utilizado: `ASM6978v1` (GCA_000069785.1 / GCF_000069785.1) — genoma completo de *C. taiwanensis* LMG 19424

```bash
mkdir -p genome
cd genome

wget https://ftp.ensemblgenomes.ebi.ac.uk/pub/bacteria/current/fasta/bacteria_105_collection/cupriavidus_taiwanensis_gca_900249975/dna/Cupriavidus_taiwanensis_gca_900249975.CBM2634_.dna.toplevel.fa.gz
gunzip Cupriavidus_taiwanensis_gca_900249975.CBM2634_.dna.toplevel.fa.gz

```

## Indexação do Genoma

```bash
bwa index Cupriavidus_taiwanensis_gca_900249975.CBM2634_.dna.toplevel.fa
```

## Alinhamento de reads

```bash
bwa mem \
Cupriavidus_taiwanensis_gca_900249975.CBM2634_.dna.toplevel.fa \SRR3943818_1.fastq \SRR3943818_2.fastq > alinhamento.sam
```
> Ferramenta utilizada bwa

## Conversão e indexação do BAM

```bash
mkdir alignment
cd alignment

samtools sort -o alinhamento_ordenado.bam alinhamento.sam
samtools index alinhamento_ordenado.bam
samtools flagstat alinhamento_ordenado.bam

```
