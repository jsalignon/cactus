
* Introduction: [Cactus](/README.md), [Quick Start](/docs/1_Intro__Quick_start.md), 
* Overview: [Graph](/docs/2_Overview__Graph.md), [Outputs structure](/docs/2_Overview__Outputs_structure.md)
* Install: [Dependencies](/docs/3_Install__Dependencies.md), [Containers](/docs/3_Install__Containers.md), [Data](/docs/3_Install__Data.md), [Test_datasets](/docs/3_Install__Test_datasets.md)
* Run: [Input Data](/docs/4_Run__Input_data.md), [Input Files](/docs/4_Run__Input_files.md), [Parameters](/docs/4_Run__Parameters.md)
* Preprocessing: ATAC: [Method](/docs/5_AP__Method.md), [Figures](/docs/5_AP__Figures.md), [MultiQC](/docs/5_AP__MultiQC.md), mRNA: [Method](/docs/6_MP__Method.md), [MultiQC](/docs/6_MP__MultiQC.md)
* Differential Abundance: [ATAC](/docs/7_DA__DiffBind.md), [mRNA](/docs/7_DA__Sleuth.md), [Figures](/docs/7_DA__Figures.md), [Tables](/docs/7_DA__Tables.md)
* Splitting peak sets: [Split](/docs/8_SP__Splitting_peak_sets.md), [Venn Diagrams](/docs/8_SP__Venn_diagrams.md)
* Enrichment: [ChIP](/docs/9_Enrich__CHIP.md), [Chromatin state](/docs/9_Enrich__Chromatin_states.md), [Motifs](/docs/9_Enrich__Motifs.md), [Func. Anno.](/docs/9_Enrich__Functional_annotations.md), [Figures](/docs/9_Enrich__Figures.md), [Tables](/docs/9_Enrich__Tables.md)


[](END_OF_MENU)







# Short Summary

The pipeline’s input are raw fastq-files and provide detailed plots and tables. ATAC-Seq data is processed following Harvard’s ATAC-seq Guidelines (the original version from 2017; link: https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html). As the pipeline was originally developed to analyze ATAC-Seq data, it contains a plethora of quality controls for checking the quality of the ATAC-Seq data. Now, though, the pipeline can be run using only mRNA-Seq data. MultiQC reports are provided for both kind of data.
For ATAC-Seq, reads are merged, trimmed and aligned to the genome via bowtie2 (ref). Aligned reads are filtered (low quality, duplicates, mitochondrial), shifted (atac shift (ref)) and extended to XX bp. Narrow peaks are then called using Macs2 (ref). These are further split, filtered (blacklist, gDNA) and annotated via CHIPseeker (ref). Finally, Differentially Accessible Regions are determined using DiffBind (ref).
For mRNA seq, transcripts are quantified using Kallisto (ref), and Differential Gene Expression Analysis is carried via the sleuth R package (ref).
Differential Abundance (DA) analysis refers to both Differentially Accessibility Analysis and Differential Gene Expression Analysis. Standardized outputs are produced for all DA results. This include Volcano plots and PCA plots, and DA results tables. DA results are subsequently filtered with user defined filters (False Discovery Rate, peak annotation, Fold Change). Split and filtered DA results are saved as tables.
Split and filtered DA results are then used as input for various enrichment analysis. For both ATAC-Seq and mRNA-Seq genes (DA peaks’ closest gene and DA genes) and genomic regions (DA peaks and DA genes’ promoter) are exported. Genes are used for ontologies and pathway enrichment analysis. Genomic regions are used for chromatin state, Transcription Factor known motifs and ENCODE CHIP enrichment analysis. Additionally, self-overlap enrichment analysis of genes and peaks is performed. Results are saved as tables, and displayed via Venn Diagrams, Barplots and Heatmaps.

# Graph

![Cactus all steps](/images/cactus_all_steps.png "Cactus all steps")



