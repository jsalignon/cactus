

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
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
- **_params.barplots__df_plots_**: An R dataframe that contains parameters to be used for each of the possible enrichment categories (i.e. data types). The default parameters (see below) can be used as a template to modify the wished parameter. Here are the parameters that can be set within this data.frame:
    
    - **_padj_threshold_**: If no adjusted pvalue is above this threshold the process is stopped and no figure is made.  
    
    - **_signed_padj_**: Should enrichment and depletion be shown (T) or enrichment only (F).  
    
    - **_add_var_**: Add a variable to the plots as a small dot. Options: 'none' (nothing added; default), 'L2OR' (log2 odd ratio), 'ov_da' (overlap of DA entries with target; i.e. counts), 'padj_loglog' (pvalues in a log scale (higher values equals lower pvalues). formula: `log10(-log10(pval) + 1)`).  
    
    - **_add_number_**: Write the number count on the plots.  
    
    - **_max_terms_**: Number of terms to display.  
    
    - **_max_characters_**: The limit of target names length. Longer targt names are cut.   
    
    - **_default parameters template_**:  
```
barplots__df_plots = 'data.frame(
data_type      = c("func_anno", "CHIP" , "motifs", "chrom_states", "genes_self", "peaks_self"),
padj_threshold = c(     0.05  ,    0.05,    0.05 ,        0.05   ,      0.05   ,       0.05  ),
signed_padj    = c(     T     ,    T   ,    T    ,        T      ,      T      ,       T     ),
add_var        = c(    "none" ,  "none",  "none" ,      "none"   ,    "none"   ,     "none"  ),
add_number     = c(     F     ,    F   ,    F    ,        F      ,      T      ,       T     ),
max_terms      = c(    30     ,   30   ,   30    ,       30      ,     30      ,      30     ),
max_characters = c(    50     ,   50   ,   50    ,       50      ,     50      ,      50     )
)'
```

### Outputs
- **Barplots**: 
  - `Figures_Individual/3_Enrichment/Barplots__${EC}/${key}__barplot.pdf` 
  - `Figures_Merged/3_Enrichment/Barplots__${EC}.pdf`.

>**_Note_:** The key for this process is `${ET}__${PA}__${FC}__${TV}__${COMP}__{EC}`.

### Examples
- **Genes self**:  
<img src="/docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__genes_self__barplot.png" width="700" />  

- **Peaks self**:   
<img src="/docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__peaks_self__barplot.png" width="700" />  

- **Functional annotations Biological Processes**:   
<img src="/docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__func_anno_BP__barplot.png" width="700" />  

- **Functional annotations KEGG**:   
<img src="/docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__func_anno_KEGG__barplot.png" width="700" />  

- **Chromatin states**:   
<img src="/docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__chrom_states__barplot.png" width="700" />  

- **CHIP**:   
<img src="/docs/examples/png/ATAC__all__down__1000__hmg4_vs_ctl__CHIP__barplot.png" width="700" />  

- **Motifs**:   
<img src="/docs/examples/png/ATAC__prom__down__200__hmg4_vs_ctl__motifs__barplot.png" width="700" />  



## Figures__making_enrichment_heatmap

### Description
This process takes as input all enrichment results for comparisons of a given group (as specified in the [comparisons.tsv file](/docs/3_Inputs/Design.md#comparisons.tsv), `${GRP}` key) and that share the same keys for `${ET}` (Experiment type), `${PA}` (Peak assignment), `${TV}` (Threshold value) and `${EC}` (Enrichment category), filters the most relevant terms, and produces a heatmap.  

The heatmap shows the selected terms on the y-axis and the comparisons with fold change type (`COMP_FC`) on the x-axis with this format: `${condition_1}_${condition_2}_${FC}`.  

The order of the `COMP_FC` entries on the x-axis, and on the y-axis for the `peaks_self`, `genes_self` enrichment categories is defined by *comparison.tsv* input file as well as the *up_down_pattern* parameter that can be set up within the *params.heatmaps__df_plots* parameter (see next subsection below).  

The order of the terms of the `chrom_states` enrichment categories (chromatin states) is defined by the chromatin state group (as defined in the original publication).  

All terms are shown for the `peaks_self`, `genes_self` and `chrom_states` enrichment categories.  

For the `CHIP`, `motifs` and `func_anno` enrichment categories a function has been created to select terms of interest (see the `params.heatmaps__df_filter_terms` parameter).
Briefly, this function first remove terms with similar names. Next, it selects the top x shared terms (significant in multiple `COMP_FC`). Then, it selects the top y terms for each `COMP_FC`. After that, the terms with the lowest pvalues accross all `COMP_FC` are selected to reach the wished number of terms. Finally, hierarchical clustering (with euclidian distance) is performed to order terms by similarity. 

Cells are colored with signed and binned adjusted pvalues as described in the [previous process](/docs/6_Enrich/Figures.md#Figures__making_enrichment_barplots) and several options are available in both processes through the *heatmaps__df_plots* parameter.  


### Parameters
- **_params.heatmaps__seed_**: random seed for the selection of terms. Default: 38.
- **_params.heatmaps__df_plots_**: An R dataframe that contains parameters to be used for each of the possible enrichment categories (i.e. data types). The default parameters (see below) can be used as a template to modify the wished parameter. Here are the parameters that can be set within this data.frame:
    
    - **_padj_threshold_**: If no adjusted pvalue is above this threshold the process is stopped and no figure is made.  
    
    - **_up_down_pattern_**: The pattern of how Fold Changes are displayed. Options: "UDUD" (up, down, up, down...) or "UUDD" (up, up, ..., down, down ...).  
    
    - **_signed_padj_**: Should enrichment and depletion be shown (T) or enrichment only (F).  
    
    - **_add_var_**: Add a variable to the plots as a small dot. Options: 'none' (nothing added; default), 'L2OR' (log2 odd ratio), 'ov_da' (overlap of DA entries with target; i.e. counts), 'padj_loglog' (pvalues in a log scale (higher values equals lower pvalues). formula: `log10(-log10(pval) + 1)`).  
    
    - **_add_number_**: Write the number count on the plots.  
    
    - **_max_characters_**: The limit of target names length. Longer targt names are cut.  
    
    - **_default parameters template_**:  
```
heatmaps__df_plots = 'data.frame(
  data_type       = c("func_anno",  "CHIP", "motifs", "chrom_states", "genes_self", "peaks_self"),
  padj_threshold  = c(     0.05  ,   0.05 ,    0.05 ,        0.05   ,      0.05   ,       0.05  ),
  up_down_pattern = c(    "UUDD" ,  "UUDD",  "UUDD" ,      "UUDD"   ,    "UUDD"   ,     "UUDD"  ),
  signed_padj     = c(     T     ,     T  ,    T    ,        T      ,      T      ,       T     ),
  add_var         = c(    "none" ,  "none",  "none" ,      "none"   ,    "none"   ,     "none"  ),
  add_number      = c(     F     ,     F  ,    F    ,        F      ,      T      ,       T     ),
  max_characters  = c(    50     ,    50  ,   50    ,       50      ,     50      ,      50     )
  )'
```

- **_params.heatmaps__df_filter_terms_**: An R data.frame that contains the parameters to use to filter the `CHIP`, `motifs`, `func_anno` enrichment categories. The default parameters (see below) can be used as a template to modify the wished parameter. Here are the parameters that can be set within this data.frame:

  - **_remove_similar_**: If true (T) entries similar names will be removed. Similar names is defined as entries that are the same before the final underscore; i.e. FOXO_L1 and FOXO_L2. For each similar entry group, the lowest pvalue of each entry is computed and the top **_remove_similar_n_** entries with the lowest pvalue are kept.  
  
  - **_n_shared_**: Number of shared terms to select. A threshold is defined with the **_threshold_type_** (options: "quantile" or "fixed" (i.e. pvalues)) and the **_threshold_value_** parameters. For each term, the number of `COMP_FC` that are below the threshold is counted. Terms are sorted by this count (with ties sorted randomly) and the top *n_shared* terms are selected.  

  - **_n_unique_**: Numbers of top terms to select. `top_N` is defined as `n_unique / n_comp` (with n_comp being the number of `COMP_FC`) rounded to the lower bound. Then for each `COMP_FC`, the `top_N` terms with the lowest pvalues are selected.

  - **_n_total_**: Total number of terms to select. This number should be higher than or equal to `n_shared + n_unique`. If the former is true, then remaining slots are taken by conditions with the lowest pvalues accross all `COMP_FC` (with ties sorted randomly).

  - **_default parameters template_**: 
```
heatmaps__df_plots = 'data.frame(
                data_type        = c("func_anno",  "CHIP"   , "motifs"   ),
                remove_similar   = c(     F     ,     T     ,    T       ),
                remove_similar_n = c(     2     ,     2     ,    2       ),
                n_shared         = c(     6     ,     8     ,    8       ),
                threshold_type   = c( "fixed"   , "quantile",  "quantile"),
                threshold_value  = c(     0.05  ,     0.25  ,    0.25    ),
                n_unique         = c(    20     ,    25     ,   25       ),
                n_total          = c(    26     ,    40     ,   40       )
                )'
```

### Outputs
- `Figures_Individual/3_Enrichment/Heatmaps__${EC}/${key}__heatmap.pdf` 
- `Figures_Merged/3_Enrichment/Heatmaps__${EC}.pdf`.

>**_Note_:** The key for this process is `${ET}__${PA}__${TV}__${GRP}__{EC}`, `${GRP}` being the current group of comparisons.


## Figures__merging_pdfs

### Description
This process uses [pdftk](https://www.pdflabs.com/tools/pdftk-the-pdf-toolkit/) to merge pdf. 

>**_Note_:** Output files and path are specified in the process where they were created.
