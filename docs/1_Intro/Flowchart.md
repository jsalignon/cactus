
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



![Cactus all steps](/docs/images/cactus_all_steps.png "Cactus all steps")

 - **Preprocessing**:
   - **ATAC-Seq**: Reads are merged, trimmed and aligned to the genome with [bowtie2]. Aligned reads are filtered (low quality, duplicates, mitochondrial, small contigs) and shifted (to account for the [9 bp offset of the transposase](https://doi.org/10.1038/nmeth.2688)). Narrow peaks are then called using [MACS2]. These are further split, filtered (blacklist, gDNA) and annotated with [ChIPseeker]. 
   - **mRNA-Seq**: Transcripts are quantified using [Kallisto].
 
 - **Differential Abundance (DA) analysis**: DA refers to both Differentially Accessibility Analysis done with [DiffBind] and to Differential Gene Expression Analysis done with [Sleuth]. Standardized plots (Volcano and PCA) and tables are produced. DA results are then filtered (significance/rank threshold) and split in subsets (experiment type (ATAC, mRNA or both), peak annotation type (promoter, gene body, distal non coding, ...), fold change type (up or down)). Tables and Venn Diagrams are made for each subset. For both ATAC-Seq and mRNA-Seq: genes (DA peaks’ closest gene and DA genes) and genomic regions (DA peaks and DA genes’ promoter) are exported. 
 
 - **Enrichment analysis**: Genes are used for ontologies and pathway enrichment analysis. Genomic regions are used for chromatin state, motifs and CHIP enrichment analysis. Additionally, self-overlap enrichment analysis of genes and peaks is performed. Results are saved as tables, and displayed as Barplots and Heatmaps.


[Bowtie2]: https://www.nature.com/articles/nmeth.1923
[ChIPseeker]: https://doi.org/10.1093/bioinformatics/btv145
[Kallisto]: https://doi.org/10.1038/nbt.3519
[Sleuth]: https://doi.org/10.1038/nmeth.4324
[MACS2]: https://doi.org/10.1101/496521 
[DiffBind]: https://doi.org/10.1038/nature10730

