

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [Data](/docs/2_Install/Data.md), [Test_datasets](/docs/2_Install/Test_datasets.md)
* [Run](/docs/3_Run/3_Run.md): [Input Data](/docs/3_Run/Input_data.md), [Input Files](/docs/3_Run/Input_files.md), [Parameters](/docs/3_Run/Parameters.md)

[](END_OF_MENU)



![Cactus all steps](/docs/images/cactus_all_steps.png "Cactus all steps")

 - ATAC-Seq preprocessing: Reads are merged, trimmed and aligned to the genome via [bowtie2]. Aligned reads are filtered (low quality, duplicates, mitochondrial, small contigs) and shifted (to account for the [9 base pair offset](https://doi.org/10.1038/nmeth.2688)). Narrow peaks are then called using [MACS2]. These are further split, filtered (blacklist, gDNA) and annotated via [ChIPseeker]. Finally, Differentially Accessible Regions are determined using [DiffBind].
 
 - mRNA-Seq preprocessing: Transcripts are quantified using [Kallisto], and Differential Gene Expression Analysis is carried via the [Sleuth] R package.
 
 - Differential Abundance (DA) analysis: DA refers to both Differentially Accessibility Analysis and Differential Gene Expression Analysis. Standardized outputs are produced for all DA results. This include Volcano plots and PCA plots, and DA results tables. DA results are subsequently filtered with user defined filters (False Discovery Rate, peak annotation, Fold Change). Split and filtered DA results are saved as tables. For both ATAC-Seq and mRNA-Seq genes (DA peaks’ closest gene and DA genes) and genomic regions (DA peaks and DA genes’ promoter) are exported.
 
 - Enrichment analysis: Genes are used for ontologies and pathway enrichment analysis, and Venn diagrams. Genomic regions are used for chromatin state, Transcription Factor known motifs and ENCODE CHIP enrichment analysis. Additionally, self-overlap enrichment analysis of genes and peaks is performed. Results are saved as tables, and displayed via Venn Diagrams, Barplots and Heatmaps.

[Bowtie2]: https://www.nature.com/articles/nmeth.1923
[ChIPseeker]: https://doi.org/10.1093/bioinformatics/btv145
[Kallisto]: https://doi.org/10.1038/nbt.3519
[Sleuth]: https://doi.org/10.1038/nmeth.4324
[MACS2]: https://doi.org/10.1101/496521 
[DiffBind]: https://doi.org/10.1038/nature10730)

