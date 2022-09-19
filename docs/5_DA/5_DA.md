

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [fastq](/docs/3_Inputs/fastq.md), [tsv](/docs/3_Inputs/tsv.md), [config](/docs/3_Inputs/config.md), [yml](/docs/3_Inputs/yml.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



![](/docs/images/5_DA.png "Differential Abundance")

This section covers Differential Abundance Analysis (DAA), that is differential binding analysis for ATAC-Seq (with [DiffBind](https://doi.org/10.1038/nature10730)) and Differential Gene Expression Analysis for mRNA-Seq (with [sleuth](http://dx.doi.org/10.1038/nmeth.4324)), as well as splitting the DAA results in subsets according to specified filters.

DAA results are processed in a standardized way to produce homogeneous plots (PCA and volcano) and tables (csv and formatted excel) for the two data types.

Finally, subsets of DAA results are extracted to allow more fine-grained analysis of the datasets. For each subset, gene sets are sent for genes-self and pathway/ontology enrichment analysis and for making Venn Diagram plots, while genomic regions are sent for peaks-self, chromatin states, motifs and CHIP enrichment analysis.  

The subsets are made using various keys (Fold change type, significance score, peak assignment (for ATAC-Seq)). Finally, if both ATAC-Seq and mRNA-Seq data are present then results of both experiment can be combined by keeping gene sets or genomic regions that pass all filters/keys in both experiments.  
