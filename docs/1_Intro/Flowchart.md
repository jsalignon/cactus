

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [Data](/docs/2_Install/Data.md), [Test_datasets](/docs/2_Install/Test_datasets.md)
* [Run](/docs/3_Run/3_Run.md): [Input Data](/docs/3_Run/Input_data.md), [Input Files](/docs/3_Run/Input_files.md), [Parameters](/docs/3_Run/Parameters.md)

[](END_OF_MENU)



![Cactus all steps](/docs/images/cactus_all_steps.png "Cactus all steps")


The pipeline’s input are raw fastq-files and provide detailed plots and tables. ATAC-Seq data is processed following Harvard’s ATAC-seq Guidelines (the original version from 2017; link: https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html). As the pipeline was originally developed to analyze ATAC-Seq data, it contains a plethora of quality controls for checking the quality of the ATAC-Seq data. Now, though, the pipeline can be run using only mRNA-Seq data. MultiQC reports are provided for both kind of data.
For ATAC-Seq, reads are merged, trimmed and aligned to the genome via bowtie2 (ref). Aligned reads are filtered (low quality, duplicates, mitochondrial), shifted (atac shift (ref)) and extended to XX bp. Narrow peaks are then called using Macs2 (ref). These are further split, filtered (blacklist, gDNA) and annotated via CHIPseeker (ref). Finally, Differentially Accessible Regions are determined using DiffBind (ref).
For mRNA seq, transcripts are quantified using Kallisto (ref), and Differential Gene Expression Analysis is carried via the sleuth R package (ref).
Differential Abundance (DA) analysis refers to both Differentially Accessibility Analysis and Differential Gene Expression Analysis. Standardized outputs are produced for all DA results. This include Volcano plots and PCA plots, and DA results tables. DA results are subsequently filtered with user defined filters (False Discovery Rate, peak annotation, Fold Change). Split and filtered DA results are saved as tables.
Split and filtered DA results are then used as input for various enrichment analysis. For both ATAC-Seq and mRNA-Seq genes (DA peaks’ closest gene and DA genes) and genomic regions (DA peaks and DA genes’ promoter) are exported. Genes are used for ontologies and pathway enrichment analysis. Genomic regions are used for chromatin state, Transcription Factor known motifs and ENCODE CHIP enrichment analysis. Additionally, self-overlap enrichment analysis of genes and peaks is performed. Results are saved as tables, and displayed via Venn Diagrams, Barplots and Heatmaps.



