

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# List of processes

  - [Tables__formatting_csv_tables](#Tables__formatting_csv_tables)
  - [Tables__merging_csv_tables](#Tables__merging_csv_tables)
  - [Tables__saving_excel_tables](#Tables__saving_excel_tables)
  - [Tables__merging_pdfs](#Tables__merging_pdfs)



## Tables__formatting_csv_tables

### Description
This process takes as input tables from many parts of the pipeline, filter them (optional) and format them (sorting by descending FDR and ascending L2FC or L2OR, rounding values to 2 digits after the decimal).  

>**_Note_:** No output section is specified here as the output and path for the formatted csv tables are specified in the process where they were created.

### Parameters
- **_params.v_fdr_thresholds_**: Vector of thresholds for filtering tables. For each data type, entries with FDR above this threhold will be removed.
Default: 'c( mRNA_detailed = 1, ATAC_detailed = 1,
                                res_simple = 1, res_filter = 1, func_anno = 1,
                                genes_self = 1, peaks_self = 1, 
                                chrom_states = 1, CHIP = 1, motifs = 1
                                )'.


## Tables__merging_csv_tables

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`


## Tables__saving_excel_tables

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`
