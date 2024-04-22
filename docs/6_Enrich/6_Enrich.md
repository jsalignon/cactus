
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Tutorial](/docs/1_Intro/tutorial.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



This section covers enrichment analysis for each Differential Analysis Subset ([DAS](/docs/5_DA/5_DA.md)).  

For a given DAS and term of interest (i.e. a pathway, a CHIP profile...), [enrichment](/docs/6_Enrich/Enrichment.md) is estimated by computing the overlap genes/genomic regions in the DAS versus the genes/genomic regions not in the DAS. Finally, a two-sided Fisher's Exact Test is performed to determine enrichment or depletion of the term since in most cases both can be biologically meaningful (see [here](https://academic.oup.com/bioinformatics/article/23/4/401/181853?login=true)).  

>**_Note_:** See part [References](/docs/2_Install/References.md) for details on how the external enrichment databases were preprocessed.  

For each DAS, enrichment results are used to produce standardized [figures](/docs/6_Enrich/Figures.md) (barplots and heatmaps) and [tables](/docs/6_Enrich/Tables.md) (csv and Excel). Here are some considerations:
 - Barplots are made showing the top n (*max_terms* option in the *params.barplots_params_\** parameters) most significantly enriched terms for a given DAS
 - Heatmaps are made showing groups of conditions (as defined in the [comparisons.tsv file](/docs/3_Inputs/Design.md#comparisons.tsv)) on the x-axis and selected terms on the y-axis. Diffent options can be used to change the selected terms.
 - Excel tables are formatted for easier and clearer scrolling by adding coloring, conditional coloring, filters, formatting header and adjusting column widths.  
 
Finally, merged files for both tables and figures are produced to allow for easier scrolling over all results.
