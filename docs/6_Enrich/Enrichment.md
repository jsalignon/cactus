
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# Table of contents

  - [Introduction](#Introduction)
  - [Enrichment__computing_functional_annotations_overlaps](#Enrichment__computing_functional_annotations_overlaps)
  - [Enrichment__computing_genes_self_overlaps](#Enrichment__computing_genes_self_overlaps)
  - [Enrichment__computing_peaks_overlaps](#Enrichment__computing_peaks_overlaps)
  - [Enrichment__computing_motifs_overlaps](#Enrichment__computing_motifs_overlaps)
  - [Enrichment__reformatting_motifs_results](#Enrichment__reformatting_motifs_results)
  - [Enrichment__computing_enrichment_pvalues](#Enrichment__computing_enrichment_pvalues)


# Introduction

In this section, the gene lists and genomic regions from the [splitting process](/docs/5_DA/Split.md#DA_split__splitting_differential_analysis_results_in_subsets) are overlapped with various databases.
Five standardized columns are made for each database:
 - `tgt`: the target against which the overlap is computed
 - `tot_tgt`: total number of target entries
 - `tot_da`: total number of entries in the DAS (Differential Analysis Subset)
 - `ov_da`: overlap of entries from the DAS and the target
 - `tot_nda`: total number of entries not in the DAS
 - `ov_nda`: overlap of entries not in the DAS and with entries in the target

>**_Note_:** Entries not in the DAS refers to all genes or detected regions (macs2 peaks or promoter) detected in the assay that are not present in the DAS.

These standardized columns are then used in subsequent process [to compute pvalues](#Enrichment__computing_enrichment_pvalues) and making [figures](/docs/6_Enrich/Figures.md) and [tables](/docs/6_Enrich/Tables.md#Tables__formatting_csv_tables).
The columns that are unique to a particular analysis are described in the corresponding process.

The keys of each subset are then augmented by adding the EC (Enrichment Category) variable. Thus the key becomes: `${ET}__${PA}__${FC}__${TV}__${COMP}__{EC}`.  
With, as defined in the [splitting process](/docs/5_DA/Split.md#DA_split__splitting_differential_analysis_results_in_subsets), the variables: 
 - ET: Experiment Type
 - PA: Peak Annotation
 - FC: Fold Change
 - TV: Theshold Value(s)
 - COMP: Comparison
 - EC: Enrichment Category

And EC can be any of these: 
 - func_anno_{BP,MF,CC,KEGG}: Ontologie databases GO_BP, GO_CC, GO_MF and KEGG
 - CHIP: Transcription factor CHIP-Seq profiles
 - chrom_states: Chromatin states from the specified chromatin state file
 - motifs: Transcription factors motifs sequences
 - peaks_self: Genomic regions subsets from the current experiment
 - genes_self: Gene list subsets from the current experiment.

*e.g.: key = ATAC__all__down__1000__hmg4_vs_ctl__func_anno_BP__enrich.*

>**_Note_:** Please see the [References section](/docs/2_Install/References.md) for details on how the external databases were downloaded and preprocessed, as well as details on the labels of the targets used in the figures and tables.

For all genomic regions enrichment analysis, the regions not in the DAS are used as a background for computing the significance of the overlaps. While for genes enrichment analysis and option is provided (*params.use_nda_as_bg_for_func_anno*) to either non DAS genes as a background or all genes in the database.


## Enrichment__computing_functional_annotations_overlaps

### Description
Overlap of gene lists with functional annotation databases is performed using [clusterProfiler](https://doi.org/10.1089/omi.2011.0118). These columns are added to the exported table: 
 - `tgt_id`: the id of the ontology
 - `genes_id`: the list of enriched genes collapsed with a "/".

### Parameters
- **_params.do_func_anno_enrichment_**: enable or disable this process. Default: true.
- **_params.use_nda_as_bg_for_func_anno_**: use non-differentially expressed genes as the background for differentially analysis. If FALSE, all genes in the database are used. Default: 'FALSE'.
- **_params.func_anno_databases_**: which database(s) to query for functional annotation enrichment analysis (KEEG, GO BP, GO CC or GO MF). Options: 'KEGG', 'CC', 'MF', 'BP'. Default: ['BP', 'KEGG']. 
- **_params.simplify_cutoff_**: [Similarity cutoff](https://rdrr.io/bioc/clusterProfiler/man/simplify-methods.html) to removed redundant go terms. Default: 0.8. 


## Enrichment__computing_genes_self_overlaps

### Description
In this process, all genes sets from subsets of the [splitting process](/docs/5_DA/Split.md#DA_split__splitting_differential_analysis_results_in_subsets) are overlapped with each other.


## Enrichment__computing_peaks_overlaps

### Description
This process takes as input genomic regions (bed files) from various sources and overlap them with genomic regions (bed files) of subsets from the [splitting process](/docs/5_DA/Split.md#DA_split__splitting_differential_analysis_results_in_subsets).  
The input genomic regions are:
 - CHIP
 - Chromatin states (hiHMM or ChromHMM) 
 - genomic regions of subsets from the splitting process -> for computing self overlap of genomic regions subsets within the experiment.

### Parameters
- **_params.chromatin_state_1_**: Chromatin state to use. Options are listed in the `references/${specie}/encode_chromatin_states_metadata.csv` file. Mandatory. No default.
- **_params.chip_ontology_**: CHIP ontology to use to filter the ENCODE CHIP files. Options are listed in the `references/${specie}/available_chip_ontology_groups.txt` file and details on the groups can be found in the file `references/${specie}/encode_chip_metadata.csv` file. Default: 'all'.


## Enrichment__computing_motifs_overlaps

### Description
This process uses [HOMER](https://doi.org/10.1016/j.molcel.2010.05.004) to compute the overlap of genomic regions of subsets in [CIS-BP motifs](https://doi.org/10.1016/j.cell.2014.08.009).

### Parameters
- **_params.do_motif_enrichment_**: enable or disable this process. Default: true.
- **_params.homer__nb_threads_**: number of threads used by Bowtie2. Default: 6.

### Outputs
- **Homer output folder**: `Processed_Data/3_Enrichment/${EC}/${key}/${key}__homer_results.txt`


## Enrichment__reformatting_motifs_results

### Description
Homver results tables are formatted in R to add the standardized columns necessary for computing pvalues.


## Enrichment__computing_enrichment_pvalues

### Description
This process takes all overlap processes, estimates significance and format tables.  

Hypergeometric minimum-likelihood two-sided p-values (`pval`) are obtained with a Fischer test in R. Two-sided Fischer tests are recommended for GO enrichment anlaysis since in most cases both enrichment and depletion can be biologically meaningful (see [reference](https://doi.org/10.1093/bioinformatics/btl633)).  
Log2 odd ratios (L2OR) is the log2 of the test's estimate.
Pvalues are then adjusted (`padj`) using Benjamini and Hochberg's [False Discovery Rate](https://doi.org/10.1093/bioinformatics/btl633).

The `pt_da` and `pt_nda` columns are added to indicate the percentage of overlap of the target with the Differential Analysis subset (DA) (`pt_da`) or non-DA (`pt_nda`) entries.  
A gene enrichment type column is added for functional annotation enrichment, to specify the gene database used.  
Results are sorted by adjusted pvalues (`padj`, descending order) and overlap of DA results (`ov_da`, ascending order).  

Finally, each elements of the key (`ET`, `PA`, `FC`, `TV`, `COMP`) are split in a separate column in the table as well as the target (`tgt`).

### Parameters
- **_params.motifs_test_type_**: The test to use for motif inputs. If 'Binomial' a two-sided binomial test is performed instead of the two-sided Fischer test. Options: 'binomial' or 'fischer' (any value). Default: 'binomial'.

### Outputs
- **Overlap tables**:
  - `Tables_Individual/3_Enrichment/${EC}/${key}__enrich.{csv,xlsx}`
  - `Tables_Merged/3_Enrichment/${EC}.{csv,xlsx}`,
