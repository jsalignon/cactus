

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# Table of contents

  - [Enrichment__computing_functional_annotations_overlaps](#Enrichment__computing_functional_annotations_overlaps)
  - [Enrichment__computing_functional_annotations_overlaps](#Enrichment__computing_functional_annotations_overlaps)
  - [Enrichment__computing_genes_self_overlaps](#Enrichment__computing_genes_self_overlaps)
  - [Enrichment__computing_peaks_overlaps](#Enrichment__computing_peaks_overlaps)
  - [Enrichment__computing_motifs_overlaps](#Enrichment__computing_motifs_overlaps)
  - [Enrichment__reformatting_motifs_results](#Enrichment__reformatting_motifs_results)
  - [Enrichment__computing_enrichment_pvalues](#Enrichment__computing_enrichment_pvalues)


# Introduction

In this section, the gene lists and genomic regions from the [splitting process](/docs/5_DA/Split.md#DA_split__splitting_differential_abundance_results_in_subsets) are overlapped with various databases.
Five standardized columns are made for each database:
 - `tgt`: the target against which the overlap is computed
 - `tot_da`: total number of differentially abundant results
 - `ov_da`: overlap of differentially abundant results with the target
 - `tot_nda`: total number of non differentially abundant results
 - `ov_nda`: overlap of non differentially abundant results with the target.  

These standardized columns are then used in subsequent process [to compute pvalues](#Enrichment__computing_enrichment_pvalues) and making [figures](/docs/6_Enrich/Figures.md) and [tables](/docs/6_Enrich/Tables.md#Tables__formatting_csv_tables).
The columns that are unique to a particular analysis are described in the corresponding process.

The keys of each subset are then augmented by adding the EC (Enrichment Category) variable. Thus the key becomes: `key="${${ET}__${PA}__${FC}__${TV}__${COMP}__{EC}__enrich}"`.  
With, as defined in the [splitting process](/docs/5_DA/Split.md#DA_split__splitting_differential_abundance_results_in_subsets), the variables: 
 - ET = Experiment Type
 - PA: Peak Assignment
 - FC: Fold Change
 - TV: Theshold Value(s)
 - COMP: Comparison.

And EC (Enrichment Category) that can be any of these: 
 - func_anno_{BP,MF,CC,KEGG}: Ontologie databases GO_BP, GO_CC, GO_MF and KEGG
 - CHIP: Transcription factor CHIP-Seq profiles
 - chrom_states: Chromatin states from the specified chromatin state file
 - motifs: Transcription factors motifs sequences
 - peaks_self: Genomic regions subsets from the current experiment
 - genes_self: Gene list subsets from the current experiment.

*e.g.: key = ATAC__all__down__1000__hmg4_vs_ctl__func_anno_BP__enrich.*




## Enrichment__computing_functional_annotations_overlaps

### Description
Overlap of gene lists with functional annotation databases is performed using [clusterProfiler](https://doi.org/10.1089/omi.2011.0118). In the exported table, the columns: 
 - `tgt` indicates the name of the ontology
 - `tgt_id` indicates the id of the ontology
 - `genes_id` indicates the list of enriched genes collapsed with a "/".

### Parameters
- **_params.do_gene_set_enrichment_**: enable or disable this process. Default: true.
- **_params.use_nda_as_bg_for_func_anno_**: use non-differentially expressed genes as the background for differentially analysis. If FALSE, all genes in the database are used. Default: 'FALSE'.
- **_params.func_anno_databases_**: which database(s) to query for functional annotation enrichment analysis. Options: 'KEGG', 'GO_CC', 'GO_MF', 'GO_BP'. Default: ['BP', 'KEGG']. 
- **_params.simplify_cutoff_**: [Similarity cutoff](https://rdrr.io/bioc/clusterProfiler/man/simplify-methods.html) to removed redundant go terms. Default: 0.8. 

### Outputs
- `Tables_Individual/3_Enrichment/${EC}/${key}__enrich.{csv,xlsx}`
- `Tables_Merged/3_Enrichment/${EC}.{csv,xlsx}`,



## Enrichment__computing_genes_self_overlaps

### Description
In this process, all genes sets from subsets of the [splitting process](/docs/5_DA/Split.md#DA_split__splitting_differential_abundance_results_in_subsets) are overlapped with each other.


## Enrichment__computing_peaks_overlaps

### Description
This process takes as input genomic regions (bed files) from various sources and overlap them with genomic regions (bed files) of subsets from the [splitting process](/docs/5_DA/Split.md#DA_split__splitting_differential_abundance_results_in_subsets).  
The input genomic regions are:
 - [ENCODE CHIP](/docs/2_Install/References.md)
 - [Chromatin states](/docs/2_Install/References.md) (hiHMM or ChromHMM) 
 - genomic regions of subsets from the splitting process -> for computing self overlap of genomic regions subsets within the experiment.

### Parameters
- **_params.chromatin_state_1_**: Chromatin state to use. Options are listed in the `references/${specie}/encode_chromatin_states_metadata.csv` file. No default.
- **_params.chip_ontology_**: CHIP ontology to use to filter the ENCODE CHIP files. Options are listed in the `references/${specie}/available_chip_ontology_groups.txt` file and details on the groups can be found in the file `references/${specie}/encode_chip_metadata.csv` file. Default: 'all'.

### Outputs
- `Tables_Individual/3_Enrichment/${EC}/${key}__enrich.{csv,xlsx}`
- `Tables_Merged/3_Enrichment/${EC}.{csv,xlsx}`


## Enrichment__computing_motifs_overlaps

### Description

### Parameters
- **_params.do_motif_enrichment_**: enable or disable this process. Default: true.
- **_params.homer__nb_threads_**: number of threads used by Bowtie2. Default: 6.


## Enrichment__reformatting_motifs_results

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- `Tables_Individual/3_Enrichment/${EC}/${key}__enrich.{csv,xlsx}`
- `Tables_Merged/3_Enrichment/${EC}.{csv,xlsx}`


## Enrichment__computing_enrichment_pvalues

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`
