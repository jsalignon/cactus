
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# List of processes

  - [DA_split__splitting_differential_analysis_results_in_subsets](#DA_split__splitting_differential_analysis_results_in_subsets)
  - [DA_split__plotting_venn_diagrams](#DA_split__plotting_venn_diagrams)


## DA_split__splitting_differential_analysis_results_in_subsets

<a href="url"> <img src="/docs/images/splitting_filters.png" width="400" /> </a>  
Diagrams showing all filters used to split peaks into subsets.  
Dotted black and grey arrows indicate respectively potential additional filters and not-showed filters. Abbreviations: FDR - False Discovery Rate, prom – promoter, distNC – distal non-coding.  

<a href="url"> <img src="/docs/images/splitting.png" width="650" /> </a>  
Diagrams illustrating the Experiment type filter.  

### Description
This process splits Differential Analysis results into subsets (i.e., DAS - Differential Analysis Subsets) in order to do enrichment analysis on many different angles and extract the most information out of the data.  
4 filters are used to split:
  - ET: Experiment Type. Can be either 'ATAC', 'mRNA', 'both', 'both_ATAC', or 'both_mRNA'.  
  - PA: Peak Assignment. Can be any combination of 'all', 'PA_gene', 'PA_interG', 'PA_prom', 'PA_5pUTR', 'PA_3pUTR', 'PA_exon', 'PA_intron', 'PA_downst', 'PA_distIn', 'PA_UTR', 'PA_TSS', 'PA_genPro', 'PA_distNC', 'PA_3kb', 'PA_8kb', 'PA_30kb'. See [DA_ATAC__saving_detailed_results_tables](/docs/5_DA/DA_ATAC.md#DA_ATAC__saving_detailed_results_tables) for details. 'all' disable this filters (all peaks are included).
  - FC: Fold Change. To split up and down-regulated results.
  - TV: Theshold Value(s). To split results by significance thresholds. 

> **_NOTE:_** The 'both*' entries indicates that the results pass the filters in both ATAC-Seq and mRNA-Seq. 'both' is used for gene lists (i.e. to find enriched ontologies), while 'both_ATAC' and 'both_mRNA' are used for genomic regions (i.e. to find enriched motifs/CHIP). 'both_ATAC' are ATAC-Seq peaks assigned to genes that are passing the filters in mRNA-Seq data as well. 'both_mRNA' are gene promoters of genes that pass the filters in mRNA-Seq and for which there are nearby ATAC-Seq peaks assigned to the same gene and that pass the filters.  

> **_NOTE:_** The process merges mRNA-Seq and ATAC-Seq results if *experiment_types = 'both'* otherwise it works on either of the two.  

Finally, a key is made, of the form `${ET}__${PA}__${FC}__${TV}__${COMP}`, with COMP indicating the comparison. This key is used to make: 
  - bed files that contain genomic regions (i.e. to find enriched motifs/CHIP)
  - R files that contain gene sets (i.e. to find enriched ontologies, for Venn diagrams plots).

In additions, two types of tables are produced: res_simple and res_filter. These two tables contain the same columns: the 5 key components (ET, PA, FC, TV and COMP), a peak_id column (Null for mRNA-Seq results), chromosome, gene name and id, pvalue and adjusted p-value and log2 fold changes. These two tables differ in their format: 
  - res_simple: each result is reported with the filters that it passes that are combined with "|" (i.e PA: 'all|prom'). This allows to quickly browse all results.
  - res_filter: only results passing filters are reported and each passed filter is on a different line (so 'all' and 'prom' would be on two different lines in the previous example). This file should be smaller as it exclude all the non-significant results.

### Parameters
- **_params.split__threshold_type_**: Defines if the threshold cuttoff is based on FDR (adjusted p-value) or rank. Options: 'FDR', 'rank'. Default: 'FDR'. 
- **_params.split__threshold_values_**: Groovy list defining the threshold cuttoff value(s). If *params.split__threshold_type = 'rank'* all entries ranked below this value will be kept (with entries ranked from lowest (rank = 1) to highest adjusted pvalues). If *params.split__threshold_type = 'FDR'* all entries with a -log10(adjusted p-value) below this threshold will be kept. e.g., *params.split__threshold_values = [ 1.3 ]* will keep all entries with an adjusted pvalue below 0.05 (i.e., -log10(0.05) = 1.30103). Multiple thresholds can be added but from the same type (FDR or rank). Default: [ 1.3 ].
- **_params.split__peak_assignment_**: Defines the peak assignment filters to use. See [DA_ATAC__saving_detailed_results_tables](/docs/5_DA/DA_ATAC.md#DA_ATAC__saving_detailed_results_tables) for options. Default: [ 'all' ].

### Outputs
- **Gene lists**: `Processed_Data/2_Differential_Analysis/DA_split__genes_rds/${key}__genes.rds`
- **Bed files**: `Processed_Data/2_Differential_Analysis/DA_split__bed_regions/${key}__regions.bed`
- **Res simple**: 
  - `Tables_Individual/2_Differential_Analysis/res_simple/${comparison}__res_simple.{csv,xlsx}`
  - `Tables_Merged/2_Differential_Analysis/res_simple.{csv,xlsx}`
<a href="url"> <img src="/docs/examples/xlsx_png/res_simple.png" width="800" /> </a>  

- **Res filter**: 
  - `Tables_Individual/2_Differential_Analysis/res_filter/${comparison}__res_filter.{csv,xlsx}`
  - `Tables_Merged/2_Differential_Analysis/res_filter.{csv,xlsx}`
<a href="url"> <img src="/docs/examples/xlsx_png/res_filter.png" width="800" /> </a>  

      

## DA_split__plotting_venn_diagrams

### Description
This process takes as input all gene lists made by the previous process for a given comparison and generates venn diagrams for gene lists that share these keys: PA (Peak Annotation), FC (Fold Change) and TV (Theshold Value).  
Two types of plots are made:
- proportional two ways venn diagrams:  ATAC-Seq vs mRNA-Seq with FC either up or down
- fixed-size four-ways venn diagrams: ATAC-Seq vs mRNA-Seq with FC up and down.
In these plots, mRNA-Seq data has an orange filling, ATAC-Seq data has a blue filling, up-regulated genes have a purple outside line and down-regulated genes have a green purple outside line.

### Outputs
- **Two-ways venn diagrams**: 
  - `Figures_Individual/2_Differential_Analysis/Venn_diagrams__two_ways/${key}__venn_up_or_down.pdf`
  - `Figures_Merged/2_Differential_Analysis/Venn_diagrams__two_ways.pdf`
<img src="/docs/examples/png/all__down__1000__hmg4_vs_ctl__venn_up_or_down.png" width="400" />  

- **Four-ways venn diagrams**: 
  - `Figures_Individual/2_Differential_Analysis/Venn_diagrams__four_ways/${key}__venn_up_and_down.pdf`
  - `Figures_Merged/2_Differential_Analysis/Venn_diagrams__four_ways.pdf`
<img src="/docs/examples/png/all__1000__hmg4_vs_ctl__venn_up_and_down.png" width="400" />  
