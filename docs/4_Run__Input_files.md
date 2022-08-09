
* Introduction: [Cactus](/README.md), [Quick Start](/docs/1_Intro__Quick_start.md), 
* Overview: [Graph](/docs/2_Overview__Graph.md), [Outputs structure](/docs/2_Overview__Outputs_structure.md)
* Install: [Dependencies](/docs/3_Install__Dependencies.md), [Containers](/docs/3_Install__Containers.md), [Data](/docs/3_Install__Data.md), [Test_datasets](/docs/3_Install__Test_datasets.md)
* Run: [Input Data](/docs/4_Run__Input_data.md), [Input Files](/docs/4_Run__Input_files.md), [Parameters](/docs/4_Run__Parameters.md)
* Preprocessing: ATAC: [Method](/docs/5_AP__Method.md), [Figures](/docs/5_AP__Figures.md), [MultiQC](/docs/5_AP__MultiQC.md), mRNA: [Method](/docs/6_MP__Method.md), [MultiQC](/docs/6_MP__MultiQC.md)
* Differential Abundance: [ATAC](/docs/7_DA__DiffBind.md), [mRNA](/docs/7_DA__Sleuth.md), [Figures](/docs/7_DA__Figures.md), [Tables](/docs/7_DA__Tables.md)
* Splitting peak sets: [Split](/docs/8_SP__Splitting_peak_sets.md), [Venn Diagrams](/docs/8_SP__Venn_diagrams.md)
* Enrichment: [ChIP](/docs/9_Enrich__CHIP.md), [Chromatin state](/docs/9_Enrich__Chromatin_states.md), [Motifs](/docs/9_Enrich__Motifs.md), [Func. Anno.](/docs/9_Enrich__Functional_annotations.md), [Figures](/docs/9_Enrich__Figures.md), [Tables](/docs/9_Enrich__Tables.md)


[](END_OF_MENU)


# Design

Inputs of the program are raw ATAC-Seq or mRNA-Seq fastq.gz files. If there are multiple files for the same conditions, the read do not necessarily need to be merged as this can be done within cactus. 
The configuration files are expected by Cactus are: atac_seq.config, mRNA_Seq.config, comparisons.config, regions_to_remove.config, grouped_conditions_for_plots.config, base.config. All are tab separated files. Below are instructions and examples on how to write these files [note: add an example per type of file]:
-	**atac_seq.config**: one line per ATAC-Seq sample. The first entry is the sample ID, and the second is the path were the files is stored. There can be multiple paths, for multiple files, in which cases these files will be merged by Cactus.
-	**mRNA_Seq.config**: one line per mRNA-Seq sample. Same formatting as for atac_seq.config
-	**regions_to_remove.config**: This file allows to remove all reads that map to the specified regions. This is useful in particular for experiments involving RNA interference, as this is known to induce a very strong sequencing signal for the locus that is repressed. The file is configurated this way:
-	**comparisons.config**: this file allows to determine which pairs of conditions will be compared to one another for Differential Abundance Analysis. It is formatted with one comparison per line with one entry per condition. 
-	**grouped_conditions_for_plots.config**: This file allows to define groups of comparisons to plot together on the heatmaps. The format is one line per group with the group ID as first entry and the comparisons as the remaining entries. The comparisons are named this way: comparison1_vs_comparison2.
-	**base.config** (optional): this file allows to overwrite Cactus default settings with custom settings. All parameters from the run.config file can be set changed here to determine how a given experiment is analyzed. See the parameters section below for more details.


# Config
