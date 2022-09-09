

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
This process produces a barplot showing the most significant results for the input subset.  

Target names are shortened and duplicate entries are removed. The top n most significant entries are kept. 
The figure made shows on the x-axis the size of the overlap and the term name on the y axis. Entries are sorted by adjusted pvalues (descending order) and overlap of DA results (ascending order). The x-axis title indicates the total number of entries in the subset (all DA entries), and for the genomic regions subsets (i.e. from bed files) the number of entries in the background (all NDA entries).  

Adjusted p-values are signed with positive values for enrichment and negative values for depletion. The signed and binned adjusted p-values are cut into 11 bins, with these cutting points (and their signed negative values): 0.2, 0.05, 1e-5 , 1e-20 , 1e-100, 0. On the figure, enrichments are depicted in green and deplection are in purple.  

Finally, it is possible to add additional colored point to the top of the bars that represent different values (*params.barplots__add_var*) and to add the overlap count (*params.barplots__add_number*).

### Parameters
- **_params.barplots__df_plots_**: An R dataframe that contains all parameters to be used for each of the possible enrichment categories (i.e. data types). The default template can be used to modify the wished parameter: 
"barplots__df_plots = 'data.frame(
  data_type      = c("func_anno",  "CHIP"   , "motifs"   , "chrom_states", "genes_self", "peaks_self"),
  padj_threshold = c(     0.05  ,     0.05  ,    0.05    ,        0.05   ,      0.05   ,       0.05  ),
  signed_padj    = c(     T     ,     T     ,    T       ,        T      ,      T      ,       T     )
  add_var        = c(    "none" ,   "none"  ,  "none"    ,      "none"   ,    "none"   ,     "none"  ),
  add_number     = c(     F     ,     F     ,    F       ,        F      ,      T      ,       T     )
  max_terms      = c(    30     ,    30     ,   30       ,       30      ,     30      ,      30     )
  max_characters = c(    50     ,    50     ,   50       ,       50      ,     50      ,      50     )
  )'"
    - **_padj_threshold_**: If no adjusted pvalue is above this threshold the process is stopped and no figure is made.
    - **_signed_padj_**: Should enrichment and depletion be shown (T) or enrichment only (F).
    - **_add_var_**: Add a variable to the plots as a small dot. Options: 'none' (nothing added; default), 'L2OR' (log2 odd ratio), 'ov_da' (overlap of DA entries with target; i.e. counts), 'padj_loglog' (pvalues in a log scale (higher values equals lower pvalues). formula: `log10(-log10(pval) + 1)`).
    - **_add_number_**: Write the number count on the plots.
    - **_max_terms_**: Number of terms to display.
    - **_max_characters_**: The limit of target names length. Longer targt names are cut. 

### Outputs
- **Barplots**: 
  - `Figures_Individual/3_Enrichment/Barplots__${EC}/${key}__barplot.pdf` 
  - `Figures_Merged/3_Enrichment/Barplots__${EC}.pdf`.

>**_Note_:** The key for this process is `${ET}__${PA}__${FC}__${TV}__${COMP}__{EC}`.


## Figures__making_enrichment_heatmap

### Description
This process takes as input all enrichment results for comparisons of a given group (as specified in the [comparisons.tsv file](/docs/3_Inputs/Design.md#comparisons.tsv), `${GRP}` key) and that share the same keys for `${ET}` (Experiment type), `${PA}` (Peak assignment), `${TV}` (Threshold value) and `${EC}` (Enrichment category), filters the most relevant terms, and produces a heatmap.  

The heatmap shows the selected terms on the y-axis and the comparisons with fold change on the x-axis with this format: `${condition_1}_${condition_2}_${FC}`. 
Cells are colored with signed and binned adjusted pvalues (as described in the [previous process](/docs/6_Enrich/Figures.md#Figures__making_enrichment_barplots)). 

For genes and peaks self-overlap only comparisons that have from the 

makes heatmaps for the groups.

key: 
'ET', 'PF', 'TV'

### Parameters
- **_params.barplots__padj_threshold_**: If no adjusted pvalue is above this threshold the process is stop and no figure is made. Default: 0.05.

### Outputs
- `Figures_Individual/3_Enrichment/Heatmaps__${EC}/${key}__heatmap.pdf` 
- `Figures_Merged/3_Enrichment/Heatmaps__${EC}.pdf`.

>**_Note_:** The key for this process is `${ET}__${PA}__${TV}__${GRP}__{EC}`, `${GRP}` being the current group of comparisons.


## Figures__merging_pdfs

### Description

### Parameters
- **_params.XX_**: AA Default: RR.

### Outputs
- **UU**: `EE`


