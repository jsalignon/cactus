
# Introduction

CACTUS (Chromatin Accessibility and Transcriptomics Unification Software) is an open-source pipeline designed to easily analyze and visualize gene expression and/or chromatin accessibility data. It can handle any of the four Encode species (human, M. musculus, D. melanogaster and C. elegans). 
The aim of the pipeline is to make it easy for the biologists, and bioinformaticians, to generate hypothesis of which molecular factors regulate gene expression and chromatin accessibility for the studied conditions. It achieves this by providing output that are both easy to use (merged pdf, formatted Excel tables, individual or separate tables/pdf) and to interpret (multiQC, volcano/PCA plots, standardized plots for mRNA and ATAC-Seq, grouped conditions for heatmap plots, …).
Finally, a key aspect of Cactus is the emphasis on efficiency and reproducibility. This is achieved via Singularity, a containerization software, and Nextflow, a dataflow-based pipeline tool that automate parallelism. 

# Pipeline Overview

[ Add the poster of the pipeline that I made. Maybe make it more fancy by putting it in Inkscape]

The pipeline’s input are raw fastq-files and provide detailed plots and tables. ATAC-Seq data is processed following Harvard’s ATAC-seq Guidelines (the original version from 2017; link: https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html). As the pipeline was originally developed to analyze ATAC-Seq data, it contains a plethora of quality controls for checking the quality of the ATAC-Seq data. Now, though, the pipeline can be run using only mRNA-Seq data. MultiQC reports are provided for both kind of data.
For ATAC-Seq, reads are merged, trimmed and aligned to the genome via bowtie2 (ref). Aligned reads are filtered (low quality, duplicates, mitochondrial), shifted (atac shift (ref)) and extended to XX bp. Narrow peaks are then called using Macs2 (ref). These are further split, filtered (blacklist, gDNA) and annotated via CHIPseeker (ref). Finally, Differentially Accessible Regions are determined using DiffBind (ref).
For mRNA seq, transcripts are quantified using Kallisto (ref), and Differential Gene Expression Analysis is carried via the sleuth R package (ref).
Differential Abundance (DA) analysis refers to both Differentially Accessibility Analysis and Differential Gene Expression Analysis. Standardized outputs are produced for all DA results. This include Volcano plots and PCA plots, and DA results tables. DA results are subsequently filtered with user defined filters (False Discovery Rate, peak annotation, Fold Change). Split and filtered DA results are saved as tables.
Split and filtered DA results are then used as input for various enrichment analysis. For both ATAC-Seq and mRNA-Seq genes (DA peaks’ closest gene and DA genes) and genomic regions (DA peaks and DA genes’ promoter) are exported. Genes are used for ontologies and pathway enrichment analysis. Genomic regions are used for chromatin state, Transcription Factor known motifs and ENCODE CHIP enrichment analysis. Additionally, self-overlap enrichment analysis of genes and peaks is performed. Results are saved as tables, and displayed via Venn Diagrams, Barplots and Heatmaps.

For more details see the publication (manuscript under preparation).




Citation

Manuscript under preparation.
