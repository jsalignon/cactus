
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



![](/docs/images/figure_1.png "Introduction")
***Overview of Cactus.***   
(a) Key features. Icons were adapted from [Servier Medical Art](https://smart.servier.com/) and the [Database Center for the life sciences/TogoTV](https://togotv.dbcls.jp/en/pics.html). (b) Simplified workflow. (c) Example of enrichment analysis performed for a gene showing an increase in both chromatin accessibility and gene expression upon treatment. Enrichment of internal GRs and GSs  indicates enrichment of GRs and GSs in other GRs and GSs generated by the pipeline. Black lines and blue circles represent DNA and nucleosomes, respectively. Orange lines represent mRNA molecules. (d) Sub-workflow showing the creation of DASs. Dotted arrows indicate optional additional filters. Abbreviations: DAR, differentially accessible region; DEG, differentially expressed gene; ChIP, ChIP-Seq binding sites; motifs, DNA binding motifs; FDR, false discovery rate; prom, promoter; distNC, distal non-coding region.

**CACTUS** (Chromatin ACcessibility and Transcriptomics Unification Software) is an mRNA-Seq and ATAC-Seq analysis pipeline that aims to assist researchers in formulating hypotheses about the molecular mechanisms regulating their conditions of interest. The pipeline does standard preprocessing and differential abundance analysis, followed by enrichment analysis using various large-scale external datasets, such as databases of gene ontologies, pathways, DNA binding motifs, CHIP-Seq profiles and chromatin states profiles. Currently, Cactus can analyze data from any of the four ENCODE/modENCODE species: *H. sapiens*, *M. musculus*, *D. melanogaster* and *C. elegans*. The pipeline is designed to be easy to use for people without bioinformatics skills, efficient and reproducible through the use of the workflow language Nextflow, and various tools managers (Singularity, Docker, Conda, Mamba), and flexible with many parameters available to customize the analysis. Output files are easy to view (e.g., multiQC, merged and individual pdfs and tables, formatted Excel tables) and interpret (e.g., standardized downstream analysis figures, customizable heatmaps).

This introductory section provides a quick overview of how Cactus works with, including:
 - A [quick start](/docs/1_Intro/Quick_start.md) guide to get started rapidly
 - A [flowchart](/docs/1_Intro/Flowchart.md) that details all the key steps of the analysis
 - An overview of the [output](/docs/1_Intro/Outputs_structure.md) files generated by Cactus 

**Reference**: The preprint titled "*Cactus: a user-friendly and reproducible ATAC-Seq and mRNA-Seq analysis pipeline for data preprocessing, differential analysis, and enrichment analysis*" can be found [here](https://www.biorxiv.org/content/10.1101/2023.05.11.540110v1).

**Licence**: This source code is released under the MIT license, included [here](LICENSE).
