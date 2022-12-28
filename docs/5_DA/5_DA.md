
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



This section covers Differential Analysis (DA), that is differential binding analysis for [ATAC-Seq](/docs/5_DA/DA_ATAC.md) (with [DiffBind](https://doi.org/10.1038/nature10730)) and Differential Gene Expression Analysis for [mRNA-Seq](/docs/5_DA/DA_mRNA.md) (with [sleuth](http://dx.doi.org/10.1038/nmeth.4324)), as well as [splitting](/docs/5_DA/Split.md) the DA results into Differential Analysis Subsets (DAS) according to specified filters.

DA results are processed to make standardized figures (PCA and volcano) and tables (csv and formatted excel) for the two data types.

Finally, DAS are extracted to allow more fine-grained analysis of the datasets. For each subset, gene sets are sent for genes-self (comparisons with other experiments) and pathway/ontology enrichment analysis and for making Venn Diagram plots, while genomic regions are sent for peaks-self, chromatin states, motifs and CHIP enrichment analysis. Genes-self and peak-self enrichment analysis consists in comparing comparisons in the current run for subsets that pass the same filters.

The subsets are made using various filtering keys (Experiment type, Fold change type, significance score, peak assignment (for ATAC-Seq)). If both ATAC-Seq and mRNA-Seq data are present then results of both experiment can be combined by keeping gene sets or genomic regions that pass all filters/keys in both experiments.  

Figures and details on how DA results are split into subsets are shown in the [split section](/docs/5_DA/Split.md).
