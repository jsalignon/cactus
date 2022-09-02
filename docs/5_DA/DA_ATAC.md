

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Plots](/docs/6_Enrich/Plots.md), [Reports](/docs/6_Enrich/Reports.md)

[](END_OF_MENU)


# List of processes

  - [DA_ATAC__doing_differential_abundance_analysis](#DA_ATAC__doing_differential_abundance_analysis)
  - [DA_ATAC__annotating_diffbind_peaks](#DA_ATAC__annotating_diffbind_peaks)
  - [DA_ATAC__plotting_differential_abundance_results](#DA_ATAC__plotting_differential_abundance_results)
  - [DA_ATAC__saving_detailed_results_tables](#DA_ATAC__saving_detailed_results_tables)


## DA_ATAC__doing_differential_abundance_analysis

### Description
This process takes as input final filtered peaks and (1 base pair) reads.  
[DiffBind](https://doi.org/10.1038/nature10730) is used to do Differential Binding analysis between two comparisons.  
Briefly, DiffBind estimates reads abundance at selected peaks of interest (i.e. the consensus peak set) and then use differential gene expression analysis tools (i.e. [DESeq2](https://doi.org/10.1186/s13059-014-0550-8) or [edgeR](https://doi.org/10.1093/bioinformatics/btp616)) to determine peaks that are differentially bound. 


### Parameters
- **_params.diffbind__min_overlap_**: Only include peaks in at least this many peaksets when generating consensus peakset. The default behavior of cactus is to include any peak from any replicate into the consensus peak set (i.e. th = 1). Non robust signal should anyway have low p-value and be filtered away in downstream analysis. See the [dba function](https://rdrr.io/bioc/DiffBind/man/dba.html) for details. Default: 1.
- **_params.diffbind__analysis_method_**: Option to use DESeq2 or edgeR for the analysis. See the [dba function](https://rdrr.io/bioc/DiffBind/man/dba.html) for details. Default: 'DBA_DESEQ2'.
- **_params.use_input_control_**: If an input control is used, grey list regions (region of high-signal in the input) will be by estimated by DiffBind via the [GreyListChIP package](10.18129/B9.bioc.GreyListChIP) and excluded from analysis. See the [DiffBind::dba.blacklist function](https://rdrr.io/bioc/DiffBind/man/dba.blacklist.html) for details. Default: false.
- **_params.diffbind__min_count_**: Minimum read count value. Any interval with fewer than this many overlapping reads will be set to have this count. See the [dba.count function](https://rdrr.io/bioc/DiffBind/man/dba.count.html) for details. Default: 0.
- **_params.diffbind__normalize_**: Normalization method to use. See the [dba.normalize function](https://rdrr.io/bioc/DiffBind/man/dba.normalize.html) for options. Default: 'DBA_NORM_RLE'.

### Outputs
- **Consensus peaks**: `${comparison}__diffbind_peaks_gr.bed` in `Processed_Data/2_Differential_Abundance/ATAC__all_peaks__bed`.
- **Diffbind object**: `${comparison}__diffbind_peaks_dbo.rds` in `Processed_Data/2_Differential_Abundance/ATAC__all_peaks__DiffBind`.
- **Read counts by replicate (GRange object)**: `${comparison}__all_peaks.rds` in `Processed_Data/2_Differential_Abundance/ATAC__all_peaks__gRange`.


## DA_ATAC__annotating_diffbind_peaks

### Description
Peaks are annotated with [ChIPseeker](http://dx.doi.org/10.1093/bioinformatics/btv145). Each peak is assigned to its closest gene using the [annotatePeak function](https://github.com/YuLab-SMU/ChIPseeker/blob/master/R/annotatePeak.R).

### Parameters
- **_params.promoter_up_diffbind_peaks_**: promoter start; upstream from TSS site.
- **_params.promoter_down_diffbind_peaks_**: promoter end; downstream from TSS site.

### Outputs
- **Annotated peaks (data.frame object)**: `${comparison}__diffb_anno_peaks_df.rds` in `Processed_Data/2_Differential_Abundance/ATAC__all_peaks__dataframe`.
- **Annotated peaks (ChIPseeker object)**: `${comparison}__diffb_anno_peaks_cs.rds` in `Processed_Data/2_Differential_Abundance/ATAC__all_peaks__ChIPseeker`.


## DA_ATAC__plotting_differential_abundance_results

### Description
This process makes standardized (i.e. similar types of plots are produced for mRNA-Seq data) PCA and volcano plots, and some other plots produced directly by DiffBind. 

### Parameters
- **_params.diffbind_plots__fdr_threshold_**: Peaks with FDR less than or equal to this value are colored in red in the volcano plot. Default: 0.05.
- **_params.diffbind_plots__top_n_labels_**: The top n peaks with lowest FDR will have their annotated gene displayed. Default: 15.

### Outputs
- **Volcano plots**: 
  - `${comparison}__ATAC_volcano.pdf` in `Figures_Individual/2_Differential_Abundance/ATAC__volcano`
  - `ATAC__volcano.pdf` in `Figures_Merged/2_Differential_Abundance`.
- **Volcano plots**: 
  - `Figures_Individual/2_Differential_Abundance/ATAC__volcano/${comparison}__ATAC_volcano.pdf`
  - `Figures_Merged/2_Differential_Abundance/ATAC__volcano.pdf`.
- **PCA plots (PC 1 and 2)**: `${comparison}__ATAC_PCA_1_2.pdf` in `Figures_Individual/2_Differential_Abundance/ATAC__PCA_1_2`.
  - top left panel: percentage of variance explained by the top 5 first principal components
  - top right panel: PCA plot for principal components 1 and 2
  - bottom panels: genes annotated to peaks that contribute the most to principal components 1 (left) and 2 (right). Color code: red or -1 indicates that the peak is a positive contributor. Blue or +1 indicates that the peak is a negative contributor. 
- **PCA plots (PC 3 and 4)**: `${comparison}__ATAC_PCA_3_4.pdf` in `Figures_Individual/2_Differential_Abundance/ATAC__PCA_3_4`.
  - same as above but for principal componenets 3 and 4
- **Other plots**: `${comparison}__ATAC_other_plots.pdf` in `Figures_Individual/2_Differential_Abundance/ATAC__other_plots`.
  - [MA plot](https://rdrr.io/bioc/DiffBind/man/dba.plotMA.html): MA and scatter plots of differential binding analysis results; using normalization factors.
  - [Heatmap plot](https://rdrr.io/bioc/DiffBind/man/dba.plotHeatmap.html): Binding site heatmap.
  - [Venn diagram](https://rdrr.io/bioc/DiffBind/man/dba.plotVenn.html): 4-way Venn diagrams showing the first 2 replicates per condition.
- **Peaks without annotations**: `${comparison}__ATAC_non_annotated_peaks.pdf` in `Processed_Data/2_Differential_Abundance/ATAC__non_annotated_peaks`.
  - Should not be many, but if there are this file can help to inspect these peaks.


## DA_ATAC__saving_detailed_results_tables

### Description
The detailed ATAC-seq results table is created in R, which includes peak name, locations, significance, log2 fold change, annotated gene name and id, annotated region, distance to TSS (Transcription Start Site), and raw counts. In addition, the following filtering columns are added: 
  - FC_{up,down}: up or down-regulated
  - PF_{3,8}kb: absolute distance of less than 3kb (kilo bases) or 8kb from the TSS
  - PF_2u1d: between 2kb upstream and 1kb downstream the TSS
  - PF_TSS: overalp with the TSS
  - PF_genProm: peak is in a genic region or in a promoter
  - PF_genic: peak is in a genic region
  - PF_prom: peak is in a promoter
  - PF_distNC: peak is in a distal intergenic region or (in an intron but not in any of these regions: promoter, 5' UTR, 3' UTR and exon). distNC stands for distal noncoding. These regions have been shown in [Daugherty *et al.*](https://doi.org/10.1101/gr.226233.117) (First ATAC-Seq paper in *C. elegans*) to be enriched in active and repressed enhancers.
These columns can all be used in the cactus configuration files to filter for peaks matching certain annotation pattern. New filtering columns could be added in the future if needed.

### Outputs
- **Table**: `${comparison}__res_detailed_atac.{csv,xlsx}` in `Tables_Individual/2_Differential_Abundance/ATAC_detailed` and (`ATAC_detailed.{csv,xlsx}` in `Tables_Merged/2_Differential_Abundance`.

