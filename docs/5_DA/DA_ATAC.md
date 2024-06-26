
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Tutorial](/docs/1_Intro/tutorial.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# List of processes

  - [DA_ATAC__doing_differential_analysis](#DA_ATAC__doing_differential_analysis)
  - [DA_ATAC__annotating_diffbind_peaks](#DA_ATAC__annotating_diffbind_peaks)
  - [DA_ATAC__plotting_differential_analysis_results](#DA_ATAC__plotting_differential_analysis_results)
  - [DA_ATAC__saving_detailed_results_tables](#DA_ATAC__saving_detailed_results_tables)


## DA_ATAC__doing_differential_analysis

### Description
This process takes as input final filtered peaks and (1 base pair) reads.  
[DiffBind](https://doi.org/10.1038/nature10730) is used to do Differential Binding analysis between two comparisons.  
Briefly, DiffBind estimates reads abundance at selected peaks of interest (i.e. the consensus peak set) and then use differential gene expression analysis tools (i.e. [DESeq2](https://doi.org/10.1186/s13059-014-0550-8) or [edgeR](https://doi.org/10.1093/bioinformatics/btp616)) to determine peaks that are differentially bound. 


### Parameters
See the function links for details and possible options. Details on the choice of default values can be found [here](https://github.com/jsalignon/cactus/blob/main/main.nf#L2474-L2530). The parameters are:
- For the [dba](https://rdrr.io/bioc/DiffBind/man/dba.html) function: 
  - **_params.diffbind__analysis_method_**: Option to use DESeq2 or edgeR for the analysis. Default: 'DBA_EDGER'.
- For edgeR analysis method:
  - **_params.diffbind__edger_tagwise_**: If using *diffbind__analysis_method = 'edgeR'* should tag-wise dispersion estimates be computed or not. See [here](https://rdrr.io/bioc/DiffBind/src/R/DBA.R) and [here](https://www.rdocumentation.org/packages/edgeR/versions/3.14.0/topics/estimateTagwiseDisp) for details. Default: 'TRUE'.
- For the [dba.blacklist](https://rdrr.io/bioc/DiffBind/man/dba.blacklist.html) function:
  - **_params.use_input_control_**: If an input control is used, grey list regions (region of high-signal in the input) will be by estimated by DiffBind via the [GreyListChIP package](10.18129/B9.bioc.GreyListChIP) and excluded from analysis. Default: false.
  - **_params.diffbind__make_grey_list_**: Should a grey list be created or not. This option can be set to 'TRUE' only if *params.use_input_control* is also *'TRUE'*. If 'TRUE', a grey list region will be created from the input control to hide hotspot regions. Default: 'FALSE'.
- For the [dba.count](https://rdrr.io/bioc/DiffBind/man/dba.count.html) function:
  - **_params.diffbind__min_overlap_**: Only include peaks in at least this many peaksets when generating consensus peakset. The default behavior of cactus is to include any peak from any replicate into the consensus peak set (i.e. th = 1). Non robust signal should anyway have low p-value and be filtered away in downstream analysis. Default: 1.
  - **_params.diffbind__score_**: Score to use in the binding affinity matrix. Raw read counts are used for analysis. This parameter only influence the counts shown in the detailled_ATAC results tables (for each individual replicates). Default: 'DBA_SCORE_NORMALIZED'.
  - **_params.diffbind__sub_control_**: Option to determine if the input control reads should be substracted to each site in each sample. Default: 'FALSE'.
  - **_params.diffbind__scale_control_**: Option to determine if reads should be scaled by library size when using the *params.diffbind__sub_control_* option. Default: 'TRUE'.
  - **_params.diffbind__min_count_**: Minimum read count value. Any interval with fewer than this many overlapping reads will be set to have this count. Default: 0.
  - **_params.diffbind__summits_**: Option to control the summit heights and locations calculated for each peak. Default: 75.
  - **_params.diffbind__filter_**: Intervals with values lower than this are excluded from analysis. Default: 1.
- For the [dba.normalize](https://rdrr.io/bioc/DiffBind/man/dba.normalize.html) function:
  - **_params.diffbind__normalization_**: Normalization method to use. Default: 'DBA_NORM_DEFAULT'.
  - **_params.diffbind__library_size_**: Method used to calculate library size. Default: 'DBA_LIBSIZE_BACKGROUND'.
  - **_params.diffbind__background_**: Should background bins be used for normalization. Can be 'FALSE', 'TRUE' (default bin size of 15000bp), or an integer (indicating the bin size). Default: 'TRUE'. 
- For the [dba.contrast](https://rdrr.io/bioc/DiffBind/man/dba.contrast.html) function:
  - **_params.diffbind__design_**: Should contrasts be specified with a formula or not. Default: 'TRUE'.

### Outputs
- **Consensus peaks**: `Processed_Data/2_Differential_Analysis/ATAC__all_peaks__bed/${comparison}__diffbind_peaks_gr.bed`.
- **Diffbind object**: `Processed_Data/2_Differential_Analysis/ATAC__all_peaks__DiffBind/${comparison}__diffbind_peaks_dbo.rds`.
- **Read counts by replicate (GRange object)**: `Processed_Data/2_Differential_Analysis/ATAC__all_peaks__gRange/${comparison}__all_peaks.rds`.


## DA_ATAC__annotating_diffbind_peaks

### Description
Peaks are annotated with [ChIPseeker](http://dx.doi.org/10.1093/bioinformatics/btv145). Each peak is assigned to its closest gene using the [annotatePeak function](https://github.com/YuLab-SMU/ChIPseeker/blob/master/R/annotatePeak.R).

### Parameters
Parameters of the [annotatePeak](https://rdrr.io/bioc/ChIPseeker/man/annotatePeak.html) function:
- **_params.chipseeker__promoter_up_**: promoter start; upstream from TSS site. Default: 1500.
- **_params.chipseeker__promoter_down_**: promoter end; downstream from TSS site. Default: 500.
- **_params.chipseeker__overlap_**: this parameter together with the *params.chipseeker__ignore_overlap* controls the genes to which peaks are assigned to. If *params.chipseeker__overlap* equals "all" and *params.chipseeker__ignore_overlap* equals 'FALSE' then if a peak overlaps to a genomic feature (i.e., exon, intron, 5'UTR, 3'UTR, CDS) it will be assigned to this gene. Otherwise, the peak will be assigned to the neighboring gene regardless of overlap with genomic features. Options: "all", "TSS". Default: 'all'.
- **_params.chipseeker__ignore_overlap_**: this parameter together with the *params.chipseeker__overlap* controls the genes to which peaks are assigned to. If *params.chipseeker__overlap* equals "all" and *params.chipseeker__ignore_overlap* equals 'FALSE' then if a peak overlaps to a genomic feature (i.e., exon, intron, 5'UTR, 3'UTR, CDS) it will be assigned to this gene. Otherwise, the peak will be assigned to the neighboring gene regardless of overlap with genomic features. Options: "all", "TSS". Default: 'FALSE'.
- **_params.chipseeker__annotation_priority_**: This parameter controls the order of priorities when there are overlaping features that overlap with the peak for assigning a genomic region for the "annotation" column. Default: "c('Promoter', '5UTR', '3UTR', 'Exon', 'Intron', 'Downstream', 'Intergenic')".
- **_params.chipseeker__ignore_upstream_**: If 'TRUE' only annotate gene at the 3' of the peak. Options: "FALSE", "TRUE". Default: 'FALSE'.
- **_params.chipseeker__ignore_downstream_**: If 'TRUE' only annotate gene at the 5' of the peak. Options: "FALSE", "TRUE". Default: 'FALSE'.

### Outputs
- **Annotated peaks (data.frame object)**: `Processed_Data/2_Differential_Analysis/ATAC__all_peaks__dataframe/${comparison}__diffb_anno_peaks_df.rds`.
- **Annotated peaks (ChIPseeker object)**: `Processed_Data/2_Differential_Analysis/ATAC__all_peaks__ChIPseeker/${comparison}__diffb_anno_peaks_cs.rds`.


## DA_ATAC__plotting_differential_analysis_results

### Description
This process makes standardized (i.e. similar types of plots are produced for mRNA-Seq data) PCA and volcano plots, boxplots of FDR by PA (Peak Assignment), and some other plots produced directly by DiffBind. 

### Parameters
- **_params.diffbind_plots__fdr_threshold_**: Peaks with FDR less than or equal to this value are colored in red in the volcano plot. Default: 0.05.
- **_params.diffbind_plots__top_n_labels_**: The top n peaks with lowest FDR will have their annotated gene displayed on the volcano plot. Default: 15.

### Outputs
- **Volcano plots**: 
  - `Figures_Individual/2_Differential_Analysis/ATAC__volcano/${comparison}__ATAC_volcano.pdf`
  - `Figures_Merged/2_Differential_Analysis/ATAC__volcano.pdf`.
<img src="/docs/examples/png/hmg4_vs_ctl__ATAC_volcano.png" width="400" />  

- **PCA plots (PC 1 and 2)**: 
  - `Figures_Individual/2_Differential_Analysis/ATAC__PCA_1_2/${comparison}__ATAC_PCA_1_2.pdf`.
  - `Figures_Merged/2_Differential_Analysis/ATAC__PCA_1_2.pdf`.
    - top left panel: percentage of variance explained by the top 5 first principal components
    - top right panel: PCA plot for principal components 1 and 2
    - bottom panels: genes annotated to peaks that contribute the most to principal components 1 (left) and 2 (right). Color code: red or -1 indicates that the peak is a positive contributor. Blue or +1 indicates that the peak is a negative contributor. 
<img src="/docs/examples/png/hmg4_vs_ctl__ATAC_PCA_1_2.png" width="400" />      

- **PCA plots (PC 3 and 4)**: 
  - `Figures_Individual/2_Differential_Analysis/ATAC__PCA_3_4/${comparison}__ATAC_PCA_3_4.pdf`.
  - `Figures_Merged/2_Differential_Analysis/ATAC__PCA_3_4.pdf`.
    - Same as above but for principal components 3 and 4.
<img src="/docs/examples/png/hmg4_vs_ctl__ATAC_PCA_3_4.png" width="400" />  

- **FDR by PA filters plots**: 
  - `Figures_Individual/2_Differential_Analysis/ATAC_FDR_by_PA/${comparison}__ATAC_FDR_by_PA.pdf`
  - `Figures_Merged/2_Differential_Analysis/ATAC_FDR_by_PA.pdf`.
  <img src="/docs/examples/png/hmg4_vs_ctl__ATAC_FDR_by_PA.png" width="400" />  

- **Other plots**; 
  - `Figures_Individual/2_Differential_Analysis/ATAC__other_plots/${comparison}__ATAC_other_plots.pdf`
  - `Figures_Merged/2_Differential_Analysis/ATAC__other_plots.pdf`
    - [MA plot](https://rdrr.io/bioc/DiffBind/man/dba.plotMA.html): MA and scatter plots of differential binding analysis results; using normalization factors.  
<a href="url"> <img src="/docs/examples/png/hmg4_vs_ctl__ATAC_other_plots-1.png" width="400" /> </a>  

    - [Heatmap plot](https://rdrr.io/bioc/DiffBind/man/dba.plotHeatmap.html): Binding site heatmap.  
<a href="url"> <img src="/docs/examples/png/hmg4_vs_ctl__ATAC_other_plots-2.png" width="400" /> </a>  

    - [Venn diagram](https://rdrr.io/bioc/DiffBind/man/dba.plotVenn.html): 4-way Venn diagrams showing the first 2 replicates per condition.
<a href="url"> <img src="/docs/examples/png/hmg4_vs_ctl__ATAC_other_plots-3.png" width="400" /> </a>  

    
- **Peaks without annotations**: 
  - `Processed_Data/2_Differential_Analysis/ATAC__non_annotated_peaks/${comparison}__ATAC_non_annotated_peaks.txt`.
    - Should not be many, but if there are this file can help to inspect these peaks.


## DA_ATAC__saving_detailed_results_tables

### Description
The detailed ATAC-seq results table is created in R, which includes peak name, coordinates, significance, log2 fold change, annotated gene name and id, annotated region, distance to TSS (Transcription Start Site), and raw counts.  

In addition, the following filtering columns are added: 
  - FC_{up,down}: up or down-regulated
  - PA_gene: genic regions annotated by ChIPseeker (i.e., intron + exon)
  - PA_interG: intergenic regions annotated by ChIPseeker
  - PA_prom: promoter regions annotated by ChIPseeker
  - PA_5pUTR: five prime UPR regions annotated by ChIPseeker
  - PA_3pUTR: three prime UPR regions annotated by ChIPseeker
  - PA_exon: exonic regions annotated by ChIPseeker
  - PA_intron: intronic regions annotated by ChIPseeker
  - PA_downst: downstream regions annotated by ChIPseeker (downstream of a gene by a maximal distance of 300 base pairs)
  - PA_distIn: distal intergenic regions annotated by ChIPseeker
  - PA_UTR: UTR regions (5pUTR + 3pUTR)
  - PA_TSS: overlap with the TSS (distanceToTSS = 0)
  - PA_genPro: genic region or promoter
  - PA_distNC: peak is in a distal intergenic region or (in an intron but not in any of these regions: promoter, 5' UTR, 3' UTR and exon). distNC stands for distal noncoding. These regions have been shown in [Daugherty *et al.*](https://doi.org/10.1101/gr.226233.117) (First ATAC-Seq paper in *C. elegans*) to be enriched in active and repressed enhancers.  
  - PA_lt{10,100,X}kb: absolute distance to the nearest gene TSS is less than 10, 100, or X kilobases, with X being defined by the parameter *params.custom_distance__less_than_X_b* (default 500 kb). Note that 10 kb is a [historically commonly used](https://www.biostars.org/p/111349) cutoff for annotating ChIP-Seq peaks, and that 100 kb and 500 kb correspond to cutoffs for proximal and distal enhancers in mouse as defined in [Xie et al](https://doi.org/10.1016/j.xgen.2023.100342).
  - PA_mt{10,100,Y}kb: absolute distance to the nearest gene TSS is more than 10, 100, or Y kilobases, with Y being defined by the parameter *params.custom_distance__more_than_Y_b* (default 500 kb).

These columns can all be used in the cactus configuration files to filter for peaks matching certain annotation pattern with the parameter *params.peak_assignment_for_splitting_subsets*. 

> **_NOTE:_** The PA_prom filter uses the `chipseeker__promoter_up` and `chipseeker__promoter_down` Cactus parameters to define the promoter regions. By default this is defined as -1500/+500 to the TSS. However, users can use this filter as an abritrary "customizable" PA filter. This can help for instance to filter out peaks that are too far away from the TSS. For instance, using the PA_prom filter with `chipseeker__promoter_up = 10000` and `chipseeker__promoter_down = 10000` would give the same result as PA_lt10kb.

> **_NOTE:_** New filtering columns could be added in the future if needed.

### Parameters
- **_params.custom_distance__less_than_X_b_**: Custom threshold for the PA_ltXkb filter used to select peaks below a given distance (in base pair) to the TSS of their closest gene. Default: 500000 (i.e., 500 kilobases).
- **_params.custom_distance__more_than_Y_b_**: Custom threshold for the PA_mtYkb filter used to select peaks above a given distance (in base pair) to the TSS of their closest gene. Default: 500000 (i.e., 500 kilobases).

### Outputs
- **Table**: 
  - `Tables_Individual/2_Differential_Analysis/ATAC_detailed/${comparison}__res_detailed_atac.{csv,xlsx}`
  - `Tables_Merged/2_Differential_Analysis/ATAC_detailed.{csv,xlsx}`.
<img src="/docs/examples/xlsx_png/ATAC_detailled_1.png" width="800" />  
<img src="/docs/examples/xlsx_png/ATAC_detailled_2.png" width="800" />  
<img src="/docs/examples/xlsx_png/ATAC_detailled_3.png" width="800" />  
