

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Plots](/docs/6_Enrich/Plots.md), [Reports](/docs/6_Enrich/Reports.md)

[](END_OF_MENU)


# List of processes

  - [DA_split__splitting_differential_abundance_results_in_subsets](#DA_split__splitting_differential_abundance_results_in_subsets)
  - [DA_split__plotting_venn_diagrams](#DA_split__plotting_venn_diagrams)


## DA_split__splitting_differential_abundance_results_in_subsets

### Description
This process splits Differential abundance results into subsets in order to do enrichment analysis on many different angle and extract the most information out of the data.  
5 filters are used to split:
  - ET: Experiment Type. Can be either 'ATAC', 'mRNA', 'both', 'both_ATAC', or 'both_mRNA'.  
  - PA: Peak Assignment. Can be any combination of 'all', 'PF_3kb', 'PF_8kb', 'PF_2u1d', 'PF_TSS', 'PF_genProm', 'PF_genic', 'PF_prom', 'PF_distNC'. See [DA_ATAC__saving_detailed_results_tables](/docs/5_DA/DA_ATAC.md#DA_ATAC__saving_detailed_results_tables) for details. 'all' disable this filters (all peaks are included). <!-- this is named Peak Filtering (PF) in main.nf for now; need to change that -->
  PF_{3,8}kb: absolute distance of less than 3kb (kilo bases) or 8kb from the TSS
  - FC: Fold Change. To split up and down-regulated results.
  - Theshold Value(s) (TV): To split results by significance thresholds. 

> **_NOTE:_** The 'both*' entries indicates that the results pass the filters in both ATAC-Seq and mRNA-Seq. 'both' is used for gene lists (i.e. to find enriched ontologies), while 'both_ATAC' and 'both_mRNA' are used for genomic regions (i.e. to find enriched motifs/CHIP). 'both_ATAC' are ATAC-Seq peaks assigned to genes that are passing the filters in mRNA-Seq data as well. 'both_mRNA' are gene promoters of genes that pass the filters in mRNA-Seq and for which there are nearby ATAC-Seq peaks assigned to the same gene and that pass the filters. 
> **_NOTE:_** The process merges mRNA-Seq and ATAC-Seq results if *experiment_types = 'both'* otherwise it works on either of the two.

Finally, a key is made, of the form `ET_${ET}__PF_${PF}__${FC}__${TV}__${COMP}`, with COMP indicating the comparison. This key is used to make: 
  - bed files that contain genomic regions (i.e. to find enriched motifs/CHIP)
  - R files that contain gene sets (i.e. to find enriched ontologies, for Venn diagrams plots).

In additions, two types of tables are produced: res_simple and res_filter. These two tables contain the same columns: the 5 key components (ET, PF, FC, TV and COMP), a peak_id column (Null for mRNA-Seq results), chromosome, gene name and id, pvalue and adjusted p-value and log2 fold changes. These two tables differ in their format: 
  - res_simple: each result is reported with the filters that it passes that are combined with "|" (i.e PF: 'all|prom'). This file allows to quickly browse all results.
  - res_filter: only results passing filters are reported and each filter it pass is on a different line (so 'all' and 'prom' would be on two different lines in the previous example). This file should be smaller as it exclude all the non-significant results.

### Parameters
- **_params.split__threshold_type_**: Defines if the threshold cuttoff is based on FDR (adjusted p-value) or rank. Default: 'FDR'. Options: 'FDR', 'rank'.
- **_params.split__threshold_values_**: Defines the threshold cuttoff value(s). Default: [ 1.3 ].
- **_params.split__peak_assignment_**: Defines the peak assignment filters to use. See [DA_ATAC__saving_detailed_results_tables](/docs/5_DA/DA_ATAC.md#DA_ATAC__saving_detailed_results_tables) for options. Default: [ 'all', 'prom', 'distNC' ].

### Outputs
- **Gene lists**: `Processed_Data/2_Differential_Abundance/DA_split__genes_rds/${key}__genes.rds`
- **Bed files**: `Processed_Data/2_Differential_Abundance/DA_split__bed_regions/${key}__regions.bed`
- **Res simple**: 
  - `Tables_Individual/2_Differential_Abundance/res_simple/${comparison}__res_simple.{csv,xlsx}`
  - `Tables_Merged/2_Differential_Abundance/res_simple.{csv,xlsx}`
- **Res filter**: 
  - `Tables_Individual/2_Differential_Abundance/res_filter/${comparison}__res_filter.{csv,xlsx}`
  - `Tables_Merged/2_Differential_Abundance/res_filter.{csv,xlsx}`

      

## DA_split__plotting_venn_diagrams

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`


