
* Introduction: [Cactus](/README.md), [Quick Start](/docs/1_Intro__Quick_start.md), 
* Overview: [Graph](/docs/2_Overview__Graph.md), [Outputs structure](/docs/2_Overview__Outputs_structure.md)
* Install: [Dependencies](/docs/3_Install__Dependencies.md), [Containers](/docs/3_Install__Containers.md), [Data](/docs/3_Install__Data.md), [Test_datasets](/docs/3_Install__Test_datasets.md)
* Run: [Input Data](/docs/4_Run__Input_data.md), [Input Files](/docs/4_Run__Input_files.md), [Parameters](/docs/4_Run__Parameters.md)
* Preprocessing: ATAC: [Method](/docs/5_AP__Method.md), [Figures](/docs/5_AP__Figures.md), [MultiQC](/docs/5_AP__MultiQC.md), mRNA: [Method](/docs/6_MP__Method.md), [MultiQC](/docs/6_MP__MultiQC.md)
* Differential Abundance: [ATAC](/docs/7_DA__DiffBind.md), [mRNA](/docs/7_DA__Sleuth.md), [Figures](/docs/7_DA__Figures.md), [Tables](/docs/7_DA__Tables.md)
* Splitting peak sets: [Split](/docs/8_SP__Splitting_peak_sets.md), [Venn Diagrams](/docs/8_SP__Venn_diagrams.md)
* Enrichment: [ChIP](/docs/9_Enrich__CHIP.md), [Chromatin state](/docs/9_Enrich__Chromatin_states.md), [Motifs](/docs/9_Enrich__Motifs.md), [Func. Anno.](/docs/9_Enrich__Functional_annotations.md), [Figures](/docs/9_Enrich__Figures.md), [Tables](/docs/9_Enrich__Tables.md)


[](END_OF_MENU)




# Cactus

CACTUS (Chromatin Accessibility and Transcriptomics Unification Software) is an open-source pipeline designed to easily analyze and visualize gene expression and/or chromatin accessibility data. It can handle any of the four Encode species (human, M. musculus, D. melanogaster and C. elegans). 
The aim of the pipeline is to make it easy for the biologists, and bioinformaticians, to generate hypothesis of which molecular factors regulate gene expression and chromatin accessibility for the studied conditions. It achieves this by providing output that are both easy to use (merged pdf, formatted Excel tables, individual or separate tables/pdf) and to interpret (multiQC, volcano/PCA plots, standardized plots for mRNA and ATAC-Seq, grouped conditions for heatmap plots, …).
Finally, a key aspect of Cactus is the emphasis on efficiency and reproducibility. This is achieved via Singularity, a containerization software, and Nextflow, a dataflow-based pipeline tool that automate parallelism. 


# Citation

Manuscript under preparation.

