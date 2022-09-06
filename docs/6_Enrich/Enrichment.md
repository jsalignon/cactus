

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# Introduction

In this section, the gene lists and genomic regions from the [splitting process](/docs/5_DA/Split.md#DA_split__splitting_differential_abundance_results_in_subsets) are overlapped with various databases.
Five standardized columns are made for each database:
 - `tgt`: the target against which the overlap is computed
 - `tot_da`: total number of differentially abundant results
 - `ov_da`: overlap of differentially abundant results with the target
 - `tot_nda`: total number of non differentially abundant results
 - `ov_nda`: overlap of non differentially abundant results with the target
These standardized columns are then used in subsequent process [to compute pvalues](/docs/6_Enrich/Overlap.md#Enrichment__computing_enrichment_pvalues) and making [report files](/docs/6_Enrich/Reports.md#Tables__formatting_csv_tables).


# List of processes

  - [Enrichment__computing_functional_annotations_overlaps](#Enrichment__computing_functional_annotations_overlaps)
  - [Enrichment__computing_genes_self_overlaps](#Enrichment__computing_genes_self_overlaps)
  - [Enrichment__computing_peaks_overlaps](#Enrichment__computing_peaks_overlaps)
  - [Enrichment__computing_motifs_overlaps](#Enrichment__computing_motifs_overlaps)
  - [Enrichment__reformatting_motifs_results](#Enrichment__reformatting_motifs_results)
  - [Enrichment__computing_enrichment_pvalues](#Enrichment__computing_enrichment_pvalues)


## Enrichment__computing_functional_annotations_overlaps

### Description
Overlap of gene lists with functional annotation databases is performed using [clusterProfiler](https://doi.org/10.1089/omi.2011.0118). In the exported table, the columns: 
 - `tgt` indicates the name of the ontology
 - `tgt_id` indicates the id of the ontology
 - `genes_id` indicates the list of enriched genes collapsed with a comma.


### Parameters
- **_params.do_gene_set_enrichment_**: enable or disable this process. Default: true.
- **_params.use_nda_as_bg_for_func_anno_**: use non-differentially expressed genes as the background for differentially analysis. Default: 'FALSE'.
- **_params.func_anno_databases_**: which databases to query for functional annotation enrichment analysis. Options: 'KEGG', 'GO_CC', 'GO_MF', 'GO_BP'. Default: ['BP', 'KEGG']. 

### Outputs
- **UU**: `EE`


## Enrichment__computing_genes_self_overlaps

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`


## Enrichment__computing_peaks_overlaps

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`


## Enrichment__computing_motifs_overlaps

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`


## Enrichment__reformatting_motifs_results

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`


## Enrichment__computing_enrichment_pvalues

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`
