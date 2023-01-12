
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# Menu

- [Configuration files](#Configuration-files)
- [Mandatory parameters](#Mandatory-parameters)
- [Global configuration file](#Global-configuration-file)
- [Design](#Design)
- [Ressources](#Ressources)
- [Mandatory parameters](#Mandatory-parameters)
- [Output Files](#Output-Files)
- [Cache](#Cache)
- [Experiment](#Experiment)
- [References](#References)
- [Processes](#Processes)
	- [1. Preprocessing: ATAC_peaks](#1-Preprocessing-ATAC_peaks)
	- [1. Preprocessing: ATAC_reads](#1-Preprocessing-ATAC_reads)
	- [1. Preprocessing: mRNA](#1-Preprocessing-mRNA)
	- [2. Differential Abundance: DA_ATAC](#2-Differential-Abundance-DA_ATAC)
	- [2. Differential Abundance: DA_mRNA](#2-Differential-Abundance-DA_mRNA)
	- [2. Differential Abundance: Split](#2-Differential-Abundance-Split)
	- [3. Enrichment: Enrichment](#3-Enrichment-Enrichment)
	- [3. Enrichment: Figures](#3-Enrichment-Figures)
	- [3. Enrichment: Tables](#3-Enrichment-Tables)
	

# Configuration files

Parameters can be set up in two different configuration files:  

- a global configuration file: that applies to all runs and with name and path: *~/.cactus.config*. An example can be found [here](/conf/.example.cactus.config).

- a run-specific configuration file: this file is the only input needed for a cactus run. It can for instance be named like that (relative path from the run folder): *yml/run.yml*. An example can be found [here](/test_datasets/worm/parameters/run.yml).

Parameters can also be set up directly on the command line. See [here](https://www.nextflow.io/docs/latest/config.html) for more details on how parameters from various sources are handled by Nextflow.


# Mandatory parameters

- **_params.species_**: species under study. Options: 'worm', 'fly', 'mouse', 'human'. Mandatory. No default.
- **_params.references_dir_**: Directory where references have been downloaded. Mandatory. No default.
- **_params.chromatin_state_**: Chromatin state to use. Options are listed in the `${params.references_dir}/${params.species}/encode_chromatin_states_metadata.csv` file. Mandatory. No default.

Additionally, these parameters are mandatory if conda, mamba or singularity is used:
- **_params.singularity_cache_dir_**: Directory where singularity images are downloaded to when Cactus is run for the first time. Mandatory if the singularity profile is used. No default.
- **_params.conda_cache_dir_**:  Directory where conda environments are downloaded to when Cactus is run for the first time. Mandatory if the conda profile is used. No default.
- **_params.mamba_cache_dir_**:  Directory where mamba environments are downloaded to when Cactus is run for the first time. Mandatory if the mamba profile is used. No default.


# Global configuration file and mandatory parameters

Any parameter can be set in the *~/.cactus.config* file.  

It is highly recommended to set up here the [mandatory parameters](/docs/3_Inputs/Parameters.md#Mandatory-parameters). One exception is `params.chromatin_state` that can be set-up globally (in *~/.cactus.config*) or locally (in the run-specific *.yml* file or on the command line) depending on the users' need.

In addition, it is recommended to set up a [NextFlow Tower token](https://www.nextflow.io/docs/latest/config.html#scope-tower) in the *~/.cactus.config* file for monitoring pipelines' execution using [Nextflow Tower](https://cloud.tower.nf/) with these parameters:
- **_params.tower_token_**: Tower token to monitor the pipeline on Tower. Default: ''. <!-- default set in conf/reports.config -->
- **_params.enable_tower_**: Directory where containers have been / will be downloaded. Default: false. <!-- default set in conf/reports.config -->
<!-- - **cactus_version**: which version of cactus to use. Default: *latest*. => To implement later!  -->


# Design

- **_params.experiment_types_**: Analyze only ATAC-Seq data, only mRNA-Seq data or both data type. Options: 'both', 'atac', 'mrna'. Default: 'both'. 
- **_params.design__mrna_fastq_**: path to the [mRNA fastq design file](/docs/3_Inputs/Design.md#ATAC-fastq). Default: 'design/mrna_fastq.tsv'.
- **_params.design__atac_fastq_**: path to the [ATAC fastq design file](/docs/3_Inputs/Design.md#mRNA-fastq). Default: 'design/atac_fastq.tsv'.
- **_params.design__comparisons_**: path to the [comparisons design file](/docs/3_Inputs/Design.md#Comparisons). Default: 'design/comparisons.tsv'.
- **_params.design__regions_to_remove_**: path to the [regions to remove design file](/docs/3_Inputs/Design.md#Regions-to-remove). Default: 'design/regions_to_remove.tsv'.
- **_params.design__genes_to_remove_**: path to the [genes to remove design file](/docs/3_Inputs/Design.md#Genes-to-remove). Default: 'design/genes_to_remove.tsv'.
- **_params.design__groups_**: path to the [groups design file](/docs/3_Inputs/Design.md#Groups). Default: 'design/groups.tsv'.
- **_params.use_input_control_**: Should a gDNA input control be used for ATAC-Seq analysis to remove [greylist regions](https://rdrr.io/bioc/DiffBind/man/dba.blacklist.html) with DiffBind, and for some quality control analysis steps. Note that the input control cannot have replicates and should have the id "input" in the *params.design__atac_fastq* file (see example [here](/test_datasets/worm/design/atac_fastq__with_input.tsv)). . Default: false. <!-- default set in conf/version.config -->


# Ressources   
<!-- default sets in conf/ressources.config -->

This part contains parameters from [Nextflow's executor scope](https://www.nextflow.io/docs/latest/config.html?highlight=queuesize#scope-executor):

- **_params.executor.queueSize_**: How many processes are queued at a given time. Default: *100*.  
- **_params.executor.$local.memory_**: Maximum total memory that will be used on the server (or local machine) during the run. Default: *80 GB*.  
- **_params.executor.$local.cpus_**: Maximum total number of CPUs that will be used on the server (or local machine) during the run. Default: *15*.  


# Output Files

- **_params.res_dir_**: Name of the directory where results will be saved. Default: 'results/Cactus_v${cactus_version}'.  <!-- default set in conf/version.config -->
- **_params.pub_mode_**: Type of publication mode to use. Options are available [here](https://www.nextflow.io/docs/latest/process.html#publishdir). Default: 'link'. <!-- default set in conf/version.config -->
- **_params.save_fastq_type_**: Saving only the last, none or all fastq files. Options: 'none', 'last', 'all'. Default: 'none'. <!-- default set in conf/run_default.config -->
- **_params.save_bam_type_**: Saving only the last, none or all bam files. Options: 'none', 'last', 'all'. Default: 'last'. <!-- default set in conf/run_default.config -->
- **_params.save_bed_type_**: Saving only the last, none or all bed files. Options: 'none', 'last', 'all'. Default: 'last'. <!-- default set in conf/run_default.config -->
- **_params.save_1bp_bam_**: Saving the 1 base pair reads after all filtering steps and tn5-shift adjustement
- adjustment the ATAC-shift. Options: 'none', 'last', 'all'. Default: false. <!-- default set in conf/run_default.config -->
- **_params.report_dir_**: Directory where reports will be saved. Default: '${params.res_dir}/Run_Info/${params.current_date}'. <!-- default set in conf/reports.config -->


# Cache 
<!-- default set in conf/version.config -->

- **_params.resume_**: Enable or disable resuming of the run with the current cache. Default: true.  
- **_params.cache_**: Type of cache to make. Options are available [here](https://www.nextflow.io/docs/latest/process.html?highlight=deep#cache). Default: 'deep'.  


# References 
<!-- default sets in conf/run_default.config -->

This part contains the path to the references. 

Cactus parses all references to simplify access to external databases to the user. However, there can be occasions where one wants to use another reference file. Any parameter from the [*species.config* file](/conf/species.config) can be modified if needed. For instance, a user analyzing worm data can try to see if human motifs are enriched by using this parameter:
```
params.pwms_motifs = "${params.references_dir}/human/homer_data/homer_motifs.txt"
```

Other species parameters that may be useful to tweak in certain situations are: *params.blacklisted_regions*, *params.encode_chip_files* or *params.chromatin_state*.


# Processes    
<!-- run this script: docs/util/get_all_parameters.sh and then use this file: docs/3_Inputs/all_config_entries.txt ; note that the last .md files need manual input for the default since they span multiple lines (for Figures.md and Tables.md) ; not also that the macs promoter parameters are duplicated 3 times and need to be removed ; also the memory_picard_ and deeptools__binsize_bigwig_creation_ parameters are duplicated and should be cleaned-->

Default parameters for the processes are defined [here](/conf/run_default.config).


## 1. Preprocessing: ATAC_peaks

- **_params.macs2__qvalue_**: q-value (minimum FDR) cutoff to call significant peaks with macs2. Default: '5e-2'.
- **_params.input_control_overlap_portion_**: sample peaks that overlap with the input control by more than this percentage (of the sample peak) will be removed. Default: 0.2.
- **_params.do_saturation_curve_**: enable or disable this process. Default: true.
Parameters of the [annotatePeak](https://rdrr.io/bioc/ChIPseeker/man/annotatePeak.html) function:
- **_params.chipseeker__promoter_up_**: promoter start; upstream from TSS site. Default: 1500.
- **_params.chipseeker__promoter_down_**: promoter end; downstream from TSS site. Default: 500.
- **_params.chipseeker__overlap_**: this parameter together with the *params.chipseeker__ignore_overlap* controls the genes to which peaks are assigned to. If *params.chipseeker__overlap* equals "all" and *params.chipseeker__ignore_overlap* equals 'FALSE' then if a peak overlaps to a genomic feature (i.e., exon, intron, 5'UTR, 3'UTR, CDS) it will be assigned to this gene. Otherwise, the peak will be assigned to the neighboring gene regardless of overlap with genomic features. Options: "all", "TSS". Default: 'all'.
- **_params.chipseeker__ignore_overlap_**: this parameter together with the *params.chipseeker__overlap* controls the genes to which peaks are assigned to. If *params.chipseeker__overlap* equals "all" and *params.chipseeker__ignore_overlap* equals 'FALSE' then if a peak overlaps to a genomic feature (i.e., exon, intron, 5'UTR, 3'UTR, CDS) it will be assigned to this gene. Otherwise, the peak will be assigned to the neighboring gene regardless of overlap with genomic features. Options: "all", "TSS". Default: 'FALSE'.
- **_params.chipseeker__annotation_priority_**: This parameter controls the order of priorities when there are overlaping features that overlap with the peak for assigning a genomic region for the "annotation" column. Default: "c('Promoter', '5UTR', '3UTR', 'Exon', 'Intron', 'Downstream', 'Intergenic')".


## 1. Preprocessing: ATAC_reads

- **_params.pigz__nb_threads_**: number of threads used for parallel compression. Default: 6.
- **_params.bowtie2__nb_threads_**: number of threads used by Bowtie2. Default: 6.
- **_params.sam_MAPQ_threshold_**: MAPQ threshold. Default: 30.
- **_params.memory_picard_**: maximum memory used by Picard. Default: '20G'.
- **_params.fastqc__nb_threads_**: number of threads used by FastQC. Default: 2.
- **_params.do_bigwig_**: enable or disable this process. Default: true.
- **_params.deeptools__binsize_bigwig_creation_**: size of the bins for the creation of the bigwig file. Smaller values increase computation time. Default: 10.
- **_params.deeptools__binsize_bigwig_correlation_**: size of the bins for computing correlation between samples. Smaller values increase computation time. Default: 10000.
- **_params.deeptools__nb_threads_**: number of threads used by DeepTools. Default: 6.
- **_params.deeptools__nb_of_1_bp_samples_**: number of 1 bp sites to sample for the coverage plots. Default: 10000.
- **_params.deeptools__normalization_method_**: normalization method to use when creating BigWig files. See [here](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html) for options. Default: 'None'.
- **_params.nb_sampled_aligned_reads_**: Number of aligned reads to sample. Default: 1000000.
- **_params.nb_sampled_trimmed_reads_**: Number of trimmed reads to sample. Default: 1000000.
- **_params.botwie2__nb_threads_**: number of threads used by Bowtie2. Default: 6.


## 1. Preprocessing: mRNA

- **_params.kallisto__nb_threads_**: number of threads used by kallisto. Default: 6.
- **_params.kallisto__bootstrap_**: Number of bootstrap samples. Default: '100'.
- **_params.kallisto__fragment_len_**: Estimated average fragment length. For single end only. Default: '180'.
- **_params.kallisto__fragment_sd_**: Estimated standard deviation of fragment length. For single end only. Default: '20'.
- **_params.fastqc__nb_threads_**: number of threads used by FastQC. Default: 2.


## 2. Differential Abundance: DA_ATAC

**Differential Binding Analysis:**  
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

**Annotations and figures:**
- Parameters of the [annotatePeak](https://rdrr.io/bioc/ChIPseeker/man/annotatePeak.html) function -> see part 1. ATAC_peaks above.
- **_params.diffbind_plots__fdr_threshold_**: Peaks with FDR less than or equal to this value are colored in red in the volcano plot. Default: 0.05.
- **_params.diffbind_plots__top_n_labels_**: The top n peaks with lowest FDR will have their annotated gene displayed on the volcano plot. Default: 15.


## 2. Differential Abundance: DA_mRNA

Figures:
- **_params.sleuth_plots__fdr_threshold_**: Peaks with FDR less than or equal to this value are colored in red in the volcano plot. Default: 0.05.
- **_params.sleuth_plots__top_n_labels_**: The top n peaks with lowest FDR will have their annotated gene displayed on the volcano plot. Default: 15.


## 2. Differential Abundance: Split

- **_params.split__threshold_type_**: Defines if the threshold cuttoff is based on FDR (adjusted p-value) or rank. Options: 'FDR', 'rank'. Default: 'FDR'. 
- **_params.split__threshold_values_**: Groovy list defining the threshold cuttoff value(s). If *params.split__threshold_type = 'rank'* all entries ranked below this value will be kept (with entries ranked from lowest (rank = 1) to highest adjusted pvalues). If *params.split__threshold_type = 'FDR'* all entries with a -log10(adjusted p-value) below this threshold will be kept. e.g., *params.split__threshold_values = [ 1.3 ]* will keep all entries with an adjusted pvalue below 0.05 (i.e., -log10(0.05) = 1.30103). Multiple thresholds can be added but from the same type (FDR or rank). Default: [ 1.3 ].
- **_params.split__peak_assignment_**: Groovy list defining the peak assignment filters to use. Options are 'all' for including all peaks, or any PA filter from the [DA_ATAC__saving_detailed_results_tables](/docs/5_DA/DA_ATAC.md#DA_ATAC__saving_detailed_results_tables) process written without the initial prefix 'PA_' (e.g., 'prom', 'distNC'...). Default: [ 'all' ].
- **_params.min_entries_DA_bed_**: Subsets with fewer entries than that will be filtered out from enrichment analysis. Default: 2. 


## 3. Enrichment: Enrichment

- **_params.disable_all_enrichments_**: If true all enrichment analysis are disabled. Default: false.
- **_params.do_only_self_enrichment_**: If true all enrichment analysis are disabled. Default: false.
- **_params.do_genes_self_enrichment_**: Enable or disable genes self enrichment analysis. Default: true.
- **_params.do_peaks_self_enrichment_**: Enable or disable peaks self enrichment analysis. Default: true.
- **_params.do_func_anno_enrichment_**: Enable or disable gene set enrichment analysis. Default: true.
- **_params.do_chrom_state_enrichment_**: Enable or disable chromatin states enrichment analysis. Default: true.
- **_params.do_chip_enrichment_**: Enable or disable CHIP-Seq enrichment analysis. Default: true.
- **_params.do_motif_enrichment_**: Enable or disable motifs enrichment analysis. Default: true.
- **_params.use_nda_as_bg_for_func_anno_**: use non-differentially expressed genes as the background for differentially analysis. If FALSE, all genes in the database are used. Default: 'FALSE'.
- **_params.func_anno_databases_**: which database(s) to query for functional annotation enrichment analysis (KEEG, GO BP, GO CC or GO MF). Options: 'KEGG', 'CC', 'MF', 'BP'. Default: ['BP', 'KEGG']. 
- **_params.simplify_cutoff_**: [Similarity cutoff](https://rdrr.io/bioc/clusterProfiler/man/simplify-methods.html) to removed redundant go terms. Default: 0.8. 
- **_params.chromatin_state_**: Chromatin state to use. Options are listed in the `${params.references_dir}/${params.species}/encode_chromatin_states_metadata.csv` file. Mandatory. No default.
- **_params.chip_ontology_**: CHIP ontology to use to filter the ENCODE CHIP files. Options are listed in the `references/${species}/available_chip_ontology_groups.txt` file and details on the groups can be found in the file `references/${species}/encode_chip_metadata.csv` file. Default: 'all'.
- **_params.homer__nb_threads_**: number of threads used by Bowtie2. Default: 6.
- **_params.motifs_test_type_**: The test to use for motif inputs. If 'Binomial' a two-sided binomial test is performed instead of the two-sided Fischer test. Options: 'binomial' or 'fischer' (any value). Default: 'binomial'.


## 3. Enrichment: Figures

- **_params.save_barplots_rds_**: Should barplots be saved as rds object or not. Default: false.
- **_params.save_heatmaps_rds_**: Should heatmaps be saved as rds object or not. Default: false.
- **_params.common_**__**_{padj_bin_breaks,barplots_params,barplots_ggplot,heatmaps_params,heatmaps_ggplot,heatmaps_filter}_**: These parameters allow to set the same parameters to each enrichment categorie. There is one parameter for each enrichment category (e.g., params.common__barplots_params). If null this parameter is disabled, otherwise the value is used as the value to set up each parameter to. Default: null.

- **_params.padj_bin_breaks_**__**_{genes_self,peaks_self,func_anno,chrom_states,CHIP,motifs}_**: A string converted to a vector in R containing the 5 adjusted p-value bins cutoff. There is one parameter for each enrichment category. Default: "c( 0.2, 0.05, 1e-5, 1e-20, 1e-100 )".
- 
- **_params.barplots_params_**__**_{genes_self,peaks_self,func_anno,chrom_states,CHIP,motifs}_**: A string converted to a vector in R containing options to customize the barplots. There is one parameter for each enrichment category. Default: "c( 0.05, T, 'none', F, 50, 30 )". The options are in order:
- **_padj_threshold_**: If no adjusted pvalue is above this threshold the process is stopped and no figure is made.
- **_signed_padj_**: Should enrichment and depletion be shown (T) or enrichment only (F).  
- **_add_var_**: Add a variable to the plots as a small dot. Options: 'none' (nothing added; default), 'L2OR' (log2 odd ratio), 'ov_da' (overlap of DA entries with target; i.e. counts), 'padj_loglog' (pvalues in a log scale (higher values equals lower pvalues). formula: `log10(-log10(pval) + 1)`).  
- **_add_number_**: Write the number count on the plots.
- **_max_characters_**: The length limit of terms names.
- **_max_terms_**: Number of terms to display.

- **_params.barplots_ggplot_**__**_{genes_self,peaks_self,func_anno,chrom_states,CHIP,motifs}_**: A string converted to a vector in R containing options to customize the appearance of the barplots by tweaking ggplot2 parameters. There is one parameter for each enrichment category. Default: "c( 11, 10, 7 )". The options are in order:
  - **_axis_text_size_**: Axis text size.
	- **_title_text_size_**: Title text size.
  - **_legend_text_size_**: Legend text size.


- **_params.heatmaps__seed_**: random seed for the selection of terms. Default: 38.

- **_params.heatmaps_params_**__**_{genes_self,peaks_self,func_anno,chrom_states,CHIP,motifs}_**: A string converted to a vector in R containing options to customize the heatmaps. There is one parameter for each enrichment category. Default for `genes_self` and `peaks_self`: "c( 0.05, T, 'none', T, 50, 'UUDD', 0 )". Default for `func_anno`, `chrom_states`, `CHIP` and `motifs`: "c( 0.05, T, 'none', F, 50, 'UUDD', 0 )". The options are in order:
  - **_padj_threshold_**: If no adjusted pvalue is above this threshold the process is stopped and no figure is made.
  - **_signed_padj_**: Should enrichment and depletion be shown (T) or enrichment only (F).
  - **_add_var_**: Add a variable to the plots as a small dot. Options: 'none' (nothing added; default), 'L2OR' (log2 odd ratio), 'ov_da' (overlap of DA entries with target; i.e. counts), 'padj_loglog' (pvalues in a log scale (higher values equals lower pvalues). formula: `log10(-log10(pval) + 1)`).
  - **_add_number_**: Write the overlap count on the cells.
  - **_max_characters_**: The limit of target names length. Longer targt names are cut.
  - **_up_down_pattern_**: The pattern of how Fold Changes are displayed. Options: "UDUD" (up, down, up, down...) or "UUDD" (up, up, ..., down, down ...).  
	- **_cell_text_size_**: Allows to control text size in the cells of the heatmap if the *add_number* parameter is set to true. If set to zero, the text size will be determined automatically by Cactus according to the number of comparisons on the heatmap.  

- **_params.heatmaps_ggplot_**__**_{genes_self,peaks_self,func_anno,chrom_states,CHIP,motifs}_**: A string converted to a vector in R containing options to customize the appearance of the heatmaps by tweaking ggplot2 parameters. There is one parameter for each enrichment category. Default: "c( 11, 10, 7 )". The options are in order:
  - **_axis_text_size_**: Axis text size.
	- **_title_text_size_**: Title text size.
  - **_legend_text_size_**: Legend text size.

- **_params.heatmaps_filter_**__**_{func_anno,CHIP,motifs}_**: A string converted to a vector in R containing options to customize the selection of terms for the heatmaps. Such filtering parameters are only available for the `func_anno`, `CHIP` and `motifs` enrichment categories. Default for `func_anno`: "c( 26, 18, 8, F, 2, 'ward.D', F)". Default for `CHIP` and `motifs`: "c( 40, 30, 10, T, 2, 'ward.D', F)". The options are in order:
  - **_n_total_**: Total number of terms to select. This number should be higher than or equal to `n_shared + n_unique`. If the former is true, then remaining slots are taken by conditions with the lowest pvalues accross all `COMP_FC` (with ties sorted randomly).
  - **_n_shared_**: Number of shared terms to select. Shared terms are defined as terms with the highest median absolute -log10 pvalue accross `COMP_FC`.
  - **_n_unique_**: Numbers of top terms to select. `top_N` is defined as `n_unique / n_comp` (with n_comp being the number of `COMP_FC`) rounded to the lower bound. Then for each `COMP_FC`, the `top_N` terms with the lowest pvalues are selected.
  - **_remove_similar_**: If true (T) entries similar names will be removed. Similar names is defined as entries that are the same before the final underscore; i.e. FOXO_L1 and FOXO_L2. For each similar entry group, the lowest pvalue of each entry is computed and the top **_remove_similar_n_** entries with the lowest pvalue are kept.
  - **_remove_similar_n_**: See *n_shared* above.
  - **_agglomeration_method_**: Agglomeration method used for hierarchical clustering of selected terms on the y-axis. See [here](https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/hclust) for options. 
  - **_select_enriched_**: Boolean indicating if only the most enriched terms should be selected (if TRUE/T) or the most enriched or depleted terms (if FALSE/F).

## 3. Enrichment: Tables

- **_params.v_fdr_thresholds_**: Vector of thresholds for filtering tables. For each data type, entries with FDR above this threhold will be removed. Default: 
```
tables__v_fdr_thresholds = 
	'c( mRNA_detailed = 1, ATAC_detailed = 1,' +
			'res_simple = 1, res_filter = 1, func_anno = 1,' +
			'genes_self = 1, peaks_self = 1, ' +
			'chrom_states = 1, CHIP = 1, motifs = 1' +
	')' 
```
- **_params.excel__add_conditional_formatting_**: To enable or disable conditional coloring. Default: 'TRUE'.
- **_params.excel__max_width_**: Maximum column width. Default: 40.
