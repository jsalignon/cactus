

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Options](/docs/3_Inputs/Options.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



![](/docs/images/6_Enrich.png "Enrichment")

This section covers enrichment analysis for each [subset from DAA](/docs/5_DA/5_DA.md).  

For a given subset and term of interest (i.e. a pathway, a CHIP profile, ...), [enrichment](/docs/6_Enrich/Enrichment.md) is measured by computing the overlap genes/genomic regions in the subset versus the genes/genomic regions not in the subset. Finally, a Fischer test is performed to determine enrichment or depletion of the term since in most cases both can be biologically meaningful (see [reference](https://academic.oup.com/bioinformatics/article/23/4/401/181853?login=true)).  

>**_Note_:** See part [References](/docs/2_Install/References.md) for details on how the CHIP, motifs and chromatin state databases were preprocessed.  

For each subset, enrichment results are used to produce standardized [figures](/docs/6_Enrich/Figures.md) (barplots and heatmaps) and [tables](/docs/6_Enrich/Tables.md) (csv and Excel). Here are some considerations:
 - barplots are made showing the top n most significantly enriched term for a given subset
 - heatmaps are made showing groups of conditions (as defined in the [comparisons.tsv file](/docs/3_Inputs/Design.md#comparisons.tsv)) on the x-axis and selected terms on the y-axis. Diffent options can be used to change the selected terms depending on the users need
 - Excel tables are formatted for easier and clearer scrolling by adding coloring, conditional coloring, filters, formatting header and adjusting column widths.  
 
Finally, merged files for both tables and figures are produced to allow for easier scrolling over all results.
