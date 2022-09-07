

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# List of processes

  - [Figures__making_enrichment_barplots](#Figures__making_enrichment_barplots)
  - [Figures__making_enrichment_heatmap](#Figures__making_enrichment_heatmap)


## Figures__making_enrichment_barplots

### Description

### Parameters
- **_params.barplots__padj_threshold_**: If no adjusted pvalue is above this threshold the process is stop and no figure is made. Default: 0.05.
- **_params.barplots__add_var_**: Add a variable to the plots as a small dot. Options: 'none' (nothing added; default), 'L2OR' (log2 odd ratio), 'ov_da' (overlap of DA entries with target; i.e. counts), 'padj_loglog' (pvalues in a log scale (higher values equals lower pvalues). formula: `log10(-log10(pval) + 1)`).
- **_params.barplots__add_number_**: Write the number count on the plots. Options: 'F' (false; default), 'T' (true).
- **_params.barplots__nb_of_terms_**: Number of terms to display. Default: 30.
- **_params.barplots__nb_of_characters_**: The limit of target names length. Longer targt names are cut. Default: 30.
30

### Outputs
- **UU**: `EE`


## Figures__making_enrichment_heatmap

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`


## Figures__merging_pdfs

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`
