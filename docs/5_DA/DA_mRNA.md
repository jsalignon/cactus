

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# List of processes

  - [DA_mRNA__doing_differential_abundance_analysis](#DA_mRNA__doing_differential_abundance_analysis)
  - [DA_mRNA__plotting_differential_abundance_results](#DA_mRNA__plotting_differential_abundance_results)
  - [DA_mRNA__saving_detailed_results_tables](#DA_mRNA__saving_detailed_results_tables)


## DA_mRNA__doing_differential_abundance_analysis

### Description
Differential gene expression analysis is carried out with [sleuth](http://dx.doi.org/10.1038/nmeth.4324).

### Outputs
- **Sleuth R object**: `1_Preprocessing/2_Differential_Abundance/mRNA__all_genes__rsleuth/${comparison}__mRNA_DEG_rsleuth.rds`
- **Promoter of all detected genes**: `1_Preprocessing/2_Differential_Abundance/mRNA__all_genes__bed_promoters/${comparison}__all_genes_prom.rds`


## DA_mRNA__plotting_differential_abundance_results

### Description
This process makes standardized (i.e. similar types of plots are produced for mRNA-Seq data) PCA and volcano plots, and some other plots produced directly by sleuth. 

### Parameters
- **_params.sleuth_plots__fdr_threshold_**: Peaks with FDR less than or equal to this value are colored in red in the volcano plot. Default: 0.05.
- **_params.sleuth_plots__top_n_labels_**: The top n peaks with lowest FDR will have their annotated gene displayed. Default: 15.

### Outputs
- **Volcano plots**: 
  - `Figures_Individual/2_Differential_Abundance/mRNA__volcano/${comparison}__mRNA_volcano.pdf`
  - `Figures_Merged/2_Differential_Abundance/mRNA__volcano.pdf`.
<img src="/docs/examples/png/hmg4_vs_ctl__mRNA_volcano.png" width="400" />      

- **PCA plots (PC 1 and 2)**: 
  - `Figures_Individual/2_Differential_Abundance/mRNA__PCA_1_2/${comparison}__mRNA_PCA_1_2.pdf`.
  - `Figures_Merged/2_Differential_Abundance/mRNA__PCA_1_2.pdf`.
    - top left panel: percentage of variance explained by the top 5 first principal components
    - top right panel: PCA plot for principal components 1 and 2
    - bottom panels: genes annotated to peaks that contribute the most to principal components 1 (left) and 2 (right). Color code: red or -1 indicates that the peak is a positive contributor. Blue or +1 indicates that the peak is a negative contributor. 
<img src="/docs/examples/png/hmg4_vs_ctl__mRNA_PCA_1_2.png" width="400" />          

- **PCA plots (PC 3 and 4)**: 
  - `Figures_Individual/2_Differential_Abundance/mRNA__PCA_3_4/${comparison}__mRNA_PCA_3_4.pdf`.
  - `Figures_Merged/2_Differential_Abundance/mRNA__PCA_3_4.pdf`.
    - Same as above but for principal components 3 and 4.
<img src="/docs/examples/png/hmg4_vs_ctl__mRNA_PCA_3_4.png" width="400" />      
    

- **Other plots**; 
  - `Figures_Individual/2_Differential_Abundance/mRNA__other_plots/${comparison}__mRNA_other_plots.pdf`
  - `Figures_Individual/2_Differential_Abundance/mRNA__other_plots/${comparison}__mRNA_other_plots.pdf`
    - [MA plot](https://rdrr.io/github/pachterlab/sleuth/man/plot_ma.html): Make an 'MA plot' for a given test. MA plots display, for each transcript, the mean of abundances across samples on the x-axis and fold change on the y-axis.  
<a href="url"> <img src="/docs/examples/png/hmg4_vs_ctl__mRNA_other_plots-1.png" width="400" /> </a>  

<img src="/docs/examples/png/hmg4_vs_ctl__mRNA_other_plots-1.png" width="400" />  
    - [Density plot](https://rdrr.io/bioc/DiffBind/man/dba.plotHeatmap.html): Plot the density of a grouping. 
<a href="url"> <img src="/docs/examples/png/hmg4_vs_ctl__mRNA_other_plots-2.png" width="400" /> </a>  


## DA_mRNA__saving_detailed_results_tables

### Description
The detailed ATAC-seq results table is created in R, which gene id, name, coordinates, significance, log2 fold change and other sleuth columns. 

### Outputs
- `Tables_Individual/2_Differential_Abundance/mRNA_detailed/${comparison}__res_detailed_mRNA.{csv,xlsx}`
- `Tables_Merged/2_Differential_Abundance/mRNA_detailed.{csv,xlsx}`.
