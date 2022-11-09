
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

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

Parameters can be passed to two different configuration files:  

- a global configuration file: that applies to all runs and with name and path: *~/.cactus.config*. An example can be found [here](/conf/.example.cactus.config).

- a run-specific configuration file: this file is the only input needed for a cactus run. It can for instance be named like that (relative path from the run folder): *yml/run.yml*. An example can be found [here](/test_datasets/worm/parameters/run.yml).


# Mandatory parameters

- **_params.references_dir_**: Directory where references have been downloaded. 
- **_params.singularity_images_dir_**: Directory where containers have been downloaded or will be downloaded (if cactus has not been run before).
- **_params.specie_**: species under study. Options: 'worm', 'fly', 'mouse', 'human'.  
- **_params.chromatin_state_**: Chromatin state to use. Options are listed in the `${params.references_dir}/${params.species}/encode_chromatin_states_metadata.csv` file.
- **_params.design__mrna_fastq_**: path to the [mRNA fastq design file](/docs/3_Inputs/Design.md#ATAC-fastq).
- **_params.design__atac_fastq_**: path to the [ATAC fastq design file](/docs/3_Inputs/Design.md#mRNA-fastq).
- **_params.design__comparisons_**: path to the [comparisons design file](/docs/3_Inputs/Design.md#Comparisons).
- **_params.design__regions_to_remove_**: path to the [regions to remove design file](/docs/3_Inputs/Design.md#Regions-to-remove).
- **_params.design__groups_**: path to the [groups design file](/docs/3_Inputs/Design.md#Groups).


# Global configuration file

Any parameter can be set in the *~/.cactus.config* file.  

However, two mandatory parameters are recommended to be put there mandatory to indicate the path where to download the references and the singularity containers.  

In addition, it is recommended to set up a [NextFlow Tower token](https://www.nextflow.io/docs/latest/config.html#scope-tower) here in order to monipor pipelines' execution using [Nextflow Tower](https://cloud.tower.nf/).  

Here are the recomended mandatory and optional parameters to put in the *~/.cactus.config* file: 
- **_params.references_dir_**: Directory where references have been downloaded. Mandatory (no default).
- **_params.singularity_images_dir_**: Directory where containers have been downloaded or will be downloaded (if cactus has not been run before). Mandatory (no default).
- **_params.tower_token_**: Tower token to monitor the pipeline on Tower. Default: ''. <!-- default set in conf/reports.config -->
- **_params.enable_tower_**: Directory where containers have been / will be downloaded. Default: false. <!-- default set in conf/reports.config -->
<!-- - **cactus_version**: which version of cactus to use. Default: *latest*. => To implement later!  -->
<!-- - **cactus_dir**: Directory where cactus is installed. Default: *~/workspace/cactus*.  => should not be needed by user... or not? -->


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
- **_params.executor.$local.cpus_**: Maximum total number of CPUs that will be used on the server (or local machine) during the run. Default: *50*.  


# Output Files

- **_params.res_dir_**: Name of the directory where results will be saved. Default: 'results/Cactus_v${cactus_version}'.  <!-- default set in conf/version.config -->
- **_params.pub_mode_**: Type of publication mode to use. Options available [here](https://www.nextflow.io/docs/latest/process.html#publishdir). Default: 'link'. <!-- default set in conf/version.config -->
- **_params.save_fastq_type_**: Saving only the last, none or all fastq files. Options: 'none', 'last', 'all'. Default: 'last'. <!-- default set in conf/run_base.config -->
- **_params.save_bam_type_**: Saving only the last, none or all bam files. Options: 'none', 'last', 'all'. Default: 'last'. <!-- default set in conf/run_base.config -->
- **_params.save_bed_type_**: Saving only the last, none or all bed files. Options: 'none', 'last', 'all'. Default: 'last'. <!-- default set in conf/run_base.config -->
- **_params.save_1bp_bam_**: Saving the 1 base pair reads after all filtering steps and tn5-shift adjustement
- adjustment the ATAC-shift. Options: 'none', 'last', 'all'. Default: 'last'. <!-- default set in conf/run_base.config -->
- **_params.report_dir_**: Directory where reports will be saved. Default: '${params.res_dir}/Run_Info/${params.current_date}'. <!-- default set in conf/reports.config -->


# Cache 
<!-- default set in conf/version.config -->

- **_params.resume_**: Enable or disable resuming of the run with the current cache. Default: true.  
- **_params.cache_**: Type of cache to make. Options avalable [here](https://www.nextflow.io/docs/latest/process.html?highlight=deep#cache). Default: 'deep'.  


# References 
<!-- default sets in conf/run_base.config -->

This part contains the path to the references. 

Cactus preparse all references to simplify access to external databases for the user. However, there can be occasions where a user want to user another reference file. Any parameter from the [*species.config* file](/conf/species.config) can be modified if needed. For instance, a user analyzing worm data can try to see if human motifs are enriched by using this parameter:
```
params.pwms_motifs = "${params.references_dir}/human/homer_data/homer_motifs.txt"
```

Other species parameters that may be useful to tweak in certain situations are: *params.blacklisted_regions*, *params.encode_chip_files* or *params.chromatin_state_1*.

The species parameter is mandatory and allows cactus to know which reference files to use:
- **_params.specie_**: species under study. Options: 'worm', 'fly', 'mouse', 'human'. Mandatory (no default).


# Processes    
<!-- run this script: docs/util/get_all_parameters.sh and then use this file: docs/3_Inputs/all_config_entries.txt ; note that the last .md files need manual input for the default since they span multiple lines (for Figures.md and Tables.md) ; not also that the macs promoter parameters are duplicated 3 times and need to be removed ; also the memory_picard_ and deeptools__binsize_bigwig_creation_ parameters are duplicated and should be cleaned-->


## 1. Preprocessing: ATAC_peaks

- **_params.macs2__qvalue_**: q-value (minimum FDR) cutoff to call significant regions. Default: '5e-2'.
- **_params.input_control_overlap_portion_**: threshold of the fraction of overlapping input control peaks to remove peaks. The percentage is regarding the treatment/sample peaks, not the input control peaks. Default: 0.2.
- **_params.do_saturation_curve_**: enable or disable this process. Default: true.
- **_params.do_raw_peak_annotation_**: to enable or disable this process. Default: true.
- **_params.macs2_peaks__promoter_up_**: promoter start; upstream from TSS site. Default: 1500.
- **_params.macs2_peaks__promoter_down_**: promoter end; downstream from TSS site. Default: 500.


## 1. Preprocessing: ATAC_reads

- **_params.pigz__nb_threads_**: number of threads used for parallel compression. Default: 6.
- **_params.bowtie2__nb_threads_**: number of threads used by Bowtie2. Default: 6.
- **_params.sam_MAPQ_threshold_**: MAPQ threshold. Default: 30.
- **_params.memory_picard_**: maximum memory used by Picard. Default: '20G'.
- **_params.fastqc__nb_threads_**: number of threads used by FastQC. Default: 2.
- **_params.do_bigwig_**: enable or disable this process. Default: true.
- **_params.deeptools__binsize_bigwig_creation_**: size of the bins in the bigwig file. Smaller values increase computation time. Default: 10000.
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

- **_params.diffbind__min_overlap_**: Only include peaks in at least this many peaksets when generating consensus peakset. The default behavior of cactus is to include any peak from any replicate into the consensus peak set (i.e. th = 1). Non robust signal should anyway have low p-value and be filtered away in downstream analysis. See the [dba function](https://rdrr.io/bioc/DiffBind/man/dba.html) for details. Default: 1.
- **_params.diffbind__analysis_method_**: Option to use DESeq2 or edgeR for the analysis. See the [dba function](https://rdrr.io/bioc/DiffBind/man/dba.html) for details. Default: 'DBA_DESEQ2'.
- **_params.use_input_control_**: If an input control is used, grey list regions (region of high-signal in the input) will be by estimated by DiffBind via the [GreyListChIP package](10.18129/B9.bioc.GreyListChIP) and excluded from analysis. See the [DiffBind::dba.blacklist function](https://rdrr.io/bioc/DiffBind/man/dba.blacklist.html) for details. Default: false.
- **_params.diffbind__min_count_**: Minimum read count value. Any interval with fewer than this many overlapping reads will be set to have this count. See the [dba.count function](https://rdrr.io/bioc/DiffBind/man/dba.count.html) for details. Default: 0.
- **_params.diffbind__normalize_**: Normalization method to use. See the [dba.normalize function](https://rdrr.io/bioc/DiffBind/man/dba.normalize.html) for options. Default: 'DBA_NORM_RLE'.
- **_params.diffbind_peaks__promoter_up_**: promoter start; upstream from TSS site. Default: 1500.
- **_params.diffbind_peaks__promoter_down_**: promoter end; downstream from TSS site. Default: 500.
- **_params.diffbind_plots__fdr_threshold_**: Peaks with FDR less than or equal to this value are colored in red in the volcano plot. Default: 0.05.
- **_params.diffbind_plots__top_n_labels_**: The top n peaks with lowest FDR will have their annotated gene displayed. Default: 15.


## 2. Differential Abundance: DA_mRNA

- **_params.sleuth_plots__fdr_threshold_**: Peaks with FDR less than or equal to this value are colored in red in the volcano plot. Default: 0.05.
- **_params.sleuth_plots__top_n_labels_**: The top n peaks with lowest FDR will have their annotated gene displayed. Default: 15.


## 2. Differential Abundance: Split

- **_params.split__threshold_type_**: Defines if the threshold cuttoff is based on FDR (adjusted p-value) or rank. Options: 'FDR', 'rank'. Default: 'FDR'. 
- **_params.split__threshold_values_**: Defines the threshold cuttoff value(s). If *params.split__threshold_type == 'rank'* all entries ranked below this value will be kept (with entries ranked from lowest (rank = 1) to highest adjusted pvalues). If *params.split__threshold_type == 'FDR'* all entries with a -log10(p-value) below this threshold will be kept. i.e. *params.split__threshold_values == [ 1.3 ]* will keep all entries with a pvalue below 0.05. Multiple thresholds can be added but from the same type (FDR or rank). Default: [ 1.3 ].
- **_params.split__peak_assignment_**: Defines the peak assignment filters to use. See [DA_ATAC__saving_detailed_results_tables](/docs/5_DA/DA_ATAC.md#DA_ATAC__saving_detailed_results_tables) for options. Default: [ 'all', 'prom', 'distNC' ].


## 3. Enrichment: Enrichment

- **_params.disable_all_enrichments_**: If true all enrichment analysis are disabled. Default: false.
- **_params.do_gene_set_enrichment_**: Enable or disable gene set enrichment analysis. Default: true.
- **_params.do_genes_self_enrichment_**: Enable or disable genes self enrichment analysis. Default: true.
- **_params.do_peaks_self_enrichment_**: Enable or disable peaks self enrichment analysis. Default: true.
- **_params.do_chrom_state_enrichment_**: Enable or disable chromatin states enrichment analysis. Default: true.
- **_params.do_motif_enrichment_**: Enable or disable motifs enrichment analysis. Default: true.
- **_params.do_chip_enrichment_**: Enable or disable CHIP-Seq enrichment analysis. Default: true.
- **_params.use_nda_as_bg_for_func_anno_**: use non-differentially expressed genes as the background for differentially analysis. If FALSE, all genes in the database are used. Default: 'FALSE'.
- **_params.func_anno_databases_**: which database(s) to query for functional annotation enrichment analysis. Options: 'KEGG', 'GO_CC', 'GO_MF', 'GO_BP'. Default: ['BP', 'KEGG']. 
- **_params.simplify_cutoff_**: [Similarity cutoff](https://rdrr.io/bioc/clusterProfiler/man/simplify-methods.html) to removed redundant go terms. Default: 0.8. 
- **_params.chromatin_state_**: Chromatin state to use. Options are listed in the `${params.references_dir}/${params.species}/encode_chromatin_states_metadata.csv` file. Mandatory (no default).
- **_params.chip_ontology_**: CHIP ontology to use to filter the ENCODE CHIP files. Options are listed in the `references/${species}/available_chip_ontology_groups.txt` file and details on the groups can be found in the file `references/${species}/encode_chip_metadata.csv` file. Default: 'all'.
- **_params.homer__nb_threads_**: number of threads used by Bowtie2. Default: 6.
- **_params.motifs_test_type_**: The test to use for motif inputs. If 'Binomial' a two-sided binomial test is performed instead of the two-sided Fischer test. Options: 'binomial' or 'fischer' (any value). Default: 'binomial'.


## 3. Enrichment: Figures

- **_params.barplots__df_plots_**: An R dataframe that contains parameters to be used for each of the possible enrichment categories (i.e. data types). The default parameters (see below) can be used as a template to modify the wished parameter. Here are the parameters that can be set within this data.frame:
    - **_padj_threshold_**: If no adjusted pvalue is above this threshold the process is stopped and no figure is made.  
    - **_signed_padj_**: Should enrichment and depletion be shown (T) or enrichment only (F).  
    - **_add_var_**: Add a variable to the plots as a small dot. Options: 'none' (nothing added; default), 'L2OR' (log2 odd ratio), 'ov_da' (overlap of DA entries with target; i.e. counts), 'padj_loglog' (pvalues in a log scale (higher values equals lower pvalues). formula: `log10(-log10(pval) + 1)`).  
    - **_add_number_**: Write the number count on the plots.  
    - **_max_terms_**: Number of terms to display.  
    - **_max_characters_**: The limit of target names length. Longer targt names are cut.   
    - **_default parameters template_**:  
```
barplots__df_plots          = 'data.frame(' +
	'data_type      = c("func_anno", "CHIP" , "motifs", "chrom_states", "genes_self", "peaks_self"),' +
	'padj_threshold = c(     0.05  ,    0.05,    0.05 ,        0.05   ,      0.05   ,       0.05  ),' +
	'signed_padj    = c(     T     ,    T   ,    T    ,        T      ,      T      ,       T     ),' +
	'add_var        = c(    "none" ,  "none",  "none" ,      "none"   ,    "none"   ,     "none"  ),' +
	'add_number     = c(     F     ,    F   ,    F    ,        F      ,      T      ,       T     ),' +
	'max_terms      = c(    30     ,   30   ,   30    ,       30      ,     30      ,      30     ),' +
	'max_characters = c(    50     ,   50   ,   50    ,       50      ,     50      ,      50     )' +
	')'
```
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
heatmaps__df_plots       = 'data.frame(' +
	'data_type       = c("func_anno",  "CHIP", "motifs", "chrom_states", "genes_self", "peaks_self"),' +
	'padj_threshold  = c(     0.05  ,   0.05 ,    0.05 ,        0.05   ,      0.05   ,       0.05  ),' +
	'up_down_pattern = c(    "UUDD" ,  "UUDD",  "UUDD" ,      "UUDD"   ,    "UUDD"   ,     "UUDD"  ),' +
	'signed_padj     = c(     T     ,     T  ,    T    ,        T      ,      T      ,       T     ),' +
	'add_var         = c(    "none" ,  "none",  "none" ,      "none"   ,    "none"   ,     "none"  ),' +
	'add_number      = c(     F     ,     F  ,    F    ,        F      ,      T      ,       T     ),' +
	'max_characters  = c(    50     ,    50  ,   50    ,       50      ,     50      ,      50     )' +
	')'
```

- **_params.heatmaps__df_filter_terms_**: An R data.frame that contains the parameters to use to filter the `CHIP`, `motifs`, `func_anno` enrichment categories. The default parameters (see below) can be used as a template to modify the wished parameter. Here are the parameters that can be set within this data.frame:
  - **_remove_similar_**: If true (T) entries similar names will be removed. Similar names is defined as entries that are the same before the final underscore; i.e. FOXO_L1 and FOXO_L2. For each similar entry group, the lowest pvalue of each entry is computed and the top **_remove_similar_n_** entries with the lowest pvalue are kept.  
  - **_n_shared_**: Number of shared terms to select. A threshold is defined with the **_threshold_type_** (options: "quantile" or "fixed" (i.e. pvalues)) and the **_threshold_value_** parameters. For each term, the number of `COMP_FC` that are below the threshold is counted. Terms are sorted by this count (with ties sorted randomly) and the top *n_shared* terms are selected.  
  - **_n_unique_**: Numbers of top terms to select. `top_N` is defined as `n_unique / n_comp` (with n_comp being the number of `COMP_FC`) rounded to the lower bound. Then for each `COMP_FC`, the `top_N` terms with the lowest pvalues are selected.
  - **_n_total_**: Total number of terms to select. This number should be higher than or equal to `n_shared + n_unique`. If the former is true, then remaining slots are taken by conditions with the lowest pvalues accross all `COMP_FC` (with ties sorted randomly).
  - **_default parameters template_**: 
```
heatmaps__df_filter_terms = 'data.frame(' +
	'data_type        = c("func_anno",  "CHIP"   , "motifs"   ),' +
	'n_shared         = c(     6     ,     8     ,    8       ),' +
	'n_unique         = c(    20     ,    25     ,   25       ),' +
	'n_total          = c(    26     ,    40     ,   40       ),' +
	'threshold_type   = c( "fixed"   , "quantile",  "quantile"),' +
	'threshold_value  = c(     0.05  ,     0.25  ,    0.25    ),' +
	'remove_similar   = c(     F     ,     T     ,    T       ),' +
	'remove_similar_n = c(     2     ,     2     ,    2       )' +
	')'
```


## 3. Enrichment: Tables

- **_params.v_fdr_thresholds_**: Vector of thresholds for filtering tables. For each data type, entries with FDR above this threhold will be removed. Default: 
```
tables__v_fdr_thresholds = 'c( mRNA_detailed = 1, ATAC_detailed = 1,' +
															'res_simple = 1, res_filter = 1, func_anno = 1,' +
															'genes_self = 1, peaks_self = 1, ' +
															'chrom_states = 1, CHIP = 1, motifs = 1' +
															')' 
```
- **_params.excel__add_conditional_formatting_**: To enable or disable conditional coloring. Default: 'TRUE'.
- **_params.excel__max_width_**: Maximum column width. Default: 40.
