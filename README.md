# Pipeline for Differential Expression (DE) analysis in R
This is the generic pipeline to perform pairwise DE analysis from the count matrix. It is designed to work with pairwise design i.e. typical Treatment Vs Control experiment type for various biological conditions. This pipeline performs following operations:
1. Read the input count matrix (pairwise design only).
2. Generate visualizations for exploratory analysis (PCA, Dendrogram, distance boxplot, sample-distance heatmap).
3. Perform DE analysis using DESeq2 and edgeR packages.
4. Determine DE genes at two significance cutoffs (FDR <= 0.01 and FDR <= 0.05) and generate output tables with annotation.
5. Additional visualizations for DE results (MA plots, Venn diagrams).
6. Generate GSEA Pre-ranked list (DESeq2 based) for GSEA analysis.

## Prerequisites:
Pipeline assume working R environment is available with following packages installed.
```
library(DESeq2)
library(ggplot2)
library(gplots)
library(tidyverse)
library(RColorBrewer)
library(pheatmap)
library(edgeR)
library(VennDiagram)
```
## Expected input:
### 1. Complete Count matrix (typically obtained from HTSeq or similar tool):

Example: The first column should have strict name **Gene_ID**, followed by control and treatment replicates denoted with numbers. 

| Gene_ID 	| control1 	| control2 	| control3 	| treated1 	| treated2 	| treated3 	|
|-----------------	|----------	|----------	|----------	|----------	|----------	|----------	|
| ENSG00000000003 	| 1013 	| 664 	| 978 	| 1090 	| 1074 	| 1025 	|
| ENSG00000000005 	| 0 	| 0 	| 0 	| 0 	| 1 	| 0 	|
| ENSG00000000419 	| 1322 	| 860 	| 1145 	| 1264 	| 1301 	| 1163 	|
| ENSG00000000457 	| 83 	| 87 	| 67 	| 94 	| 165 	| 105 	|
| ENSG00000000460 	| 305 	| 262 	| 258 	| 261 	| 376 	| 237 	|
| ENSG00000000938 	| 2 	| 0 	| 0 	| 0 	| 0 	| 1 	|


### 2. Annotation file (typically downloaded from the ENSEMBL or similar database):

Example:

| Gene_ID 	| Gene type 	| Gene name 	| Gene description 	|
|-----------------	|------------------------	|------------	|------------------------------------------------------------------------------------	|
| ENSG00000229147 	| unprocessed_pseudogene 	| SMPD4P2 	| sphingomyelin phosphodiesterase 4 pseudogene 2 [Source:HGNC Symbol;Acc:HGNC:39674] 	|
| ENSG00000256453 	| protein_coding 	| DND1 	| DND microRNA-mediated repression inhibitor 1 [Source:HGNC Symbol;Acc:HGNC:23799] 	|
| ENSG00000185813 	| protein_coding 	| PCYT2 	| phosphate cytidylyltransferase 2, ethanolamine [Source:HGNC Symbol;Acc:HGNC:8756] 	|
| ENSG00000268861 	| protein_coding 	| AC008878.3 	| Rho/Rac guanine nucleotide exchange factor 18 [Source:NCBI gene;Acc:23370] 	|
| ENSG00000281782 	| unprocessed_pseudogene 	| AC093642.6 	| F-box protein 25 (FBXO25) pseudogene 	|
| ENSG00000176749 	| protein_coding 	| CDK5R1 	| cyclin dependent kinase 5 regulatory subunit 1 [Source:HGNC Symbol;Acc:HGNC:1775] 	|

## Quick start:
Clone of download the repository:
```
git clone https://github.com/sagarutturkar/RNASeq_DE_Pipeline.git
```
**Syntax** (for Rstudio)
```
system("RScript DE_pairwise_pipeline.R   <Count_Matrix>  <Control_Name>  <Treatment_Name>  <Number_of_Control_replicates>  <Number_of_Treatment_replicates>  <Annotation.TXT>")

system("RScript Make_venn.R  DESeq2_FDR005_filtered.tsv  edgeR_FDR005_filtered.tsv DESeq2  edgeR  <custom_TAG> <Annotation.TXT>")
```

**Working Example** (for Rstudio)
```
system("RScript DE_pairwise_pipeline.R   C1.TXT  control  treated  3  3  Annotation.TXT")

system("RScript Make_venn.R  DESeq2_FDR005_filtered.tsv  edgeR_FDR005_filtered.tsv DESeq2  edgeR  FDR005   Annotation.TXT")
```

**On Linux** (command line):
```
cd test_data

Rscript ../DE_pairwise_pipeline.R C1.txt  control  treated  3  3  Annotation.txt

Rscript ../Make_Venn.R DESeq2_FDR005_filtered.tsv  edgeR_FDR005_filtered.tsv DESeq2  edgeR  FDR005   Annotation.txt
```

## Output files summary:
### Note:
> Result tables are TAB delimited TEXT files and best viewed when opened in Excel.

| File Type 	| File Name 	| File Description 	|
|---------------------------------	|-----------------------------------	|---------------------------------------------------------------------------	|
| Log Files 	| sessioninfo.txt 	| Session Info for the current run 	|
|  	| edgeR_log.txt 	| log from the edgeR package 	|
|  	| DESeq2_log.txt 	| log from the DESeq2 package 	|
| Rdata Files 	| edgeR_results.RData 	| edgeR results stored as Rdata Object 	|
|  	| DESeq2_results.RData 	| DESeq2 results stored as Rdata Object 	|
| Exploratory data visualizations 	| DESeq2_Cluster_Dendrogram.png 	| cluster dendrogram for control and treatment samples 	|
|  	| DESeq2_Cooks_distance_boxplot.png 	| heatmap showing the Euclidian distances between samples 	|
|  	| Sample_distance_matrix.png 	| Boxplot showing the distribution of Cook’s distances for each library 	|
|  	| DESeq2_PCA_with_Labels.png 	| PCA plot with lables 	|
|  	| edgeR_mds.png 	| PCA plot from edgeR 	|
| DESeq2 Result Tables 	| DESeq2_All.tsv 	| Complete Table from DESeq2 (no filter) 	|
|  	| DESeq2_FDR001_filtered.tsv 	| DESeq2 significant DE genes (FDR <= 0.01) 	|
|  	| DESeq2_FDR005_filtered.tsv 	| DESeq2 significant DE genes (FDR <= 0.05) 	|
| edgeR Result Tables 	| edgeR_All.tsv 	| Complete Table from edgeR (no filter) 	|
|  	| edgeR_FDR001_filtered.tsv 	| edgeR significant DE genes (FDR <= 0.01) 	|
|  	| edgeR_FDR005_filtered.tsv 	| edgeR significant DE genes (FDR <= 0.05) 	|
| DE Result visualizations 	| edgeR_smear.png 	| MA plot from edgeR 	|
|  	| DESeq2_MA.png 	| MA plot from DESeq2 	|
| Comaprison of DESeq2 and edgeR 	| FDR005_Overlap.tsv 	| Table for the overlapping DE genes (FDR <= 0.05) 	|
|  	| FDR005_Venn.png 	| Venn diagram for overlap between DESeq2 and edgeR DE genes at FDR <= 0.05 	|
| GSEA Pre-ranked file 	| GSEA.rnk 	| GSEA pre-ranked file 	|

## Example data:
Eample data is provided in the directory **test_data**. It contains a random counts matrix and annotations (Human).

Steps:
1. Download the two R scripts **(DE_pairwise_pipeline.R and Make_Venn.R)**.
2. Download the test_data directory.
3. Run the R scripts as below:

```
system("RScript DE_pairwise_pipeline.R   C1.TXT  control  treated  3  3  Annotation.TXT")

system("RScript Make_venn.R  DESeq2_FDR005_filtered.tsv  edgeR_FDR005_filtered.tsv DESeq2  edgeR   FDR005   Annotation.TXT")

system("RScript Make_venn.R  DESeq2_FDR001_filtered.tsv  edgeR_FDR001_filtered.tsv DESeq2  edgeR   FDR001   Annotation.TXT")
```
**Runtime:** On standrad laptop with 1 processor, the run with test data should complete under 10 minutes.
