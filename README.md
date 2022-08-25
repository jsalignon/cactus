

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [DBA](/docs/5_DA/DBA.md), [DGEA](/docs/5_DA/DGEA.md), [Split](/docs/5_DA/Split.md), [Outputs](/docs/5_DA/Outputs.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Outputs](/docs/6_Enrich/Outputs.md)

[](END_OF_MENU)






![](/docs/images/1_Intro.png "Introduction")

**_CACTUS_**: (Chromatin Accessibility and Transcriptomics Unification Software) is an open-source pipeline designed to easily analyze and visualize gene expression and/or chromatin accessibility data. It can handle any of the four Encode species (human, M. musculus, D. melanogaster and C. elegans). 
The aim of the pipeline is to make it easy for the biologists, and bioinformaticians, to generate hypothesis of which molecular factors regulate gene expression and chromatin accessibility for the studied conditions. It achieves this by providing output that are both easy to use (merged pdf, formatted Excel tables, individual or separate tables/pdf) and to interpret (multiQC, volcano/PCA plots, standardized plots for mRNA and ATAC-Seq, grouped conditions for heatmap plots, …).
Finally, a key aspect of Cactus is the emphasis on efficiency and reproducibility. This is achieved via Singularity, a containerization software, and Nextflow, a dataflow-based pipeline tool that automate parallelism. 

**_Introductory sections_**:
 - a [quick start](/docs/1_Intro/Quick_start.md) section to get started rapidly
 - a [flowchart](/docs/1_Intro/Flowchart.md) that detail all key steps of the analysis pipeline
 - an overview of the [outputs](/docs/1_Intro/Outputs_structure.md) generated by cactus 

**_Citation_**: Manuscript under preparation.

