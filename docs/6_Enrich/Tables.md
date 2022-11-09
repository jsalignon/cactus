
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
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
This process takes as input tables from many parts of the pipeline, filter them (optional) and format them (sorting by descending FDR and ascending L2FC (Log2 Fold Change) or L2OR (Log2 Odd Ratio), rounding values to 2 digits after the decimal).  

>**_Note_:** Output files and path are specified in the process where they were created.

### Parameters
- **_params.v_fdr_thresholds_**: Vector of thresholds for filtering tables. For each data type, entries with FDR above this threhold will be removed. Default: 
        'c( mRNA_detailed = 1, ATAC_detailed = 1,
            res_simple = 1, res_filter = 1, func_anno = 1,
            genes_self = 1, peaks_self = 1, 
            chrom_states = 1, CHIP = 1, motifs = 1
            )'.  


## Tables__merging_csv_tables

### Description
This process merge csv tables in R and format them (sorting by descending FDR and ascending L2FC or L2OR, rounding values to 2 digits after the decimal).

>**_Note_:** Output files and path are specified in the process where they were created.


## Tables__saving_excel_tables

### Description
This process format tables in Excel with these steps:
 - Homogeneous coloring: coloring cells by column types: with the subset keys in red, the target in green, pvalue columns in orange, L2OR and L2FC in purple, DA (Differential Abundance) columns in gold, NDA (non-DA) columns in blue, and total entry in target (tot_tgt) in grey. Headers are colored with a darker color than other cells (body) in the column.
 - Conditional coloring for the body of: 
   - adjusted pvalues: gradient from light orange (highest pvalues) to dark orange (lowest pvalues).
   - L2OR and L2FC: gradient from red (highest values) to white (zero) to blue (lowest values). Infinite values are assigned a value of +/-e99 and given the same value as the highest/lowest non-infinite value for coloring.
 - Filters are added to all columns and column width and header height are adjusted.

>**_Note_:** Output files and path are specified in the process where they were created.

### Parameters
- **_params.excel__add_conditional_formatting_**: To enable or disable conditional coloring. Default: 'TRUE'.
- **_params.excel__max_width_**: Maximum column width. Default: 40.

### Examples
- **Genes self**:  
<img src="/docs/examples/xlsx_png/genes_self.png" width="700" />  

- **Peaks self**:   
<img src="/docs/examples/xlsx_png/peaks_self.png" width="700" />  

- **Functional annotations GO-BP**:   
<img src="/docs/examples/xlsx_png/func_anno__BP.png" width="700" />  

- **Functional annotations KEGG**:   
<img src="/docs/examples/xlsx_png/func_anno__KEGG.png" width="700" />  

- **Chromatin states**:   
<img src="/docs/examples/xlsx_png/chrom_states.png" width="700" />  

- **CHIP**:   
<img src="/docs/examples/xlsx_png/CHIP.png" width="700" />  

- **Motifs**:  
<img src="/docs/examples/xlsx_png/motifs.png" width="700" />  
