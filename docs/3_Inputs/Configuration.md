

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Plots](/docs/6_Enrich/Plots.md), [Reports](/docs/6_Enrich/Reports.md)

[](END_OF_MENU)


# Configuration files

-	**_~/.cactus.config_** (mandatory): global profile that applies to all runs
	
-	**_conf/run.config_** (optional): this file allows to overwrite Cactus default settings with custom settings. All parameters from the run.config file can be set here to determine how a given experiment is analyzed. 


# Global parameters

Any parameter listed below can be set in the **_~/.cactus.config_** file. However, two parameters are mandatory to indicate the path where to download: the references (*params.references_dir*) and the singularity containers (*singularity_containers_path*). 

In addition, it is recommended to set up a NextFlow Tower token here in order to monipor pipelines' execution here, using the [tower scope](https://www.nextflow.io/docs/latest/config.html#scope-tower). 

Here is an example of a basic **_.cactus.config_** file: 

```
params.references_dir       = '/home/user/workspace/cactus/references'
singularity_containers_path = '/home/user/workspace/singularity_containers'

tower {
  accessToken = my_token
  enabled     = true
}
``` 


# List of all parameters

Below are all parameters that can be set in the **_.config_** files. The different options are listed in italic, with the first one being the default parameter.

## Tools
 - **cactus_version**: *latest*.  
Indicates which version of cactus to use.
 - **cactus_dir**: *~/workspace/cactus*.  
Directory where cactus is installed.
 - **singularity_containers_path**: *~/workspace/singularity_containers*.  
Directory where containers are downloaded.

## Ressources
- **executor.queueSize**: *50*.  
Maximum total number of CPUs that will be used on the server (or local machine) during the run.
- **executor.$local.memory**: *80 GB*s.  
Maximum total memory that will be used on the server (or local machine) during the run.
 - **nb_threads**: *6*.  
Number of CPUs to allocate by task for various steps in the analysis (kallisto, bowtie2, homer, pigz and Deeptools).
 - **memory_picard**: *20G*.  
Memory allocation for picard tools.

## Output Files
 - **out_dir**: *results/${cactus_version}*.  
Name of the directory where results will be saved.
 - **pub_mode**: *link, symlink, rellink, link, copy, copyNoFollow,move*.  
Type of [publication mode](https://www.nextflow.io/docs/latest/process.html#publishdir).
- **save_fastq_type**: *'none', 'last', 'all'*  
Saving only the last, none or all fastq files.
- **save_bam_type**: *'last', 'none', 'all'*  
Saving only the last, none or all bam files.
- **save_bed_type**: *'last', 'none', 'all'*  
Saving only the last, none or all bed files.

## Do or skip analysis
Booleans that indicates if certain parts of the pipeline should be run or not.
 - **do_motif_enrichment**: *true, false*
 - **do_chip_enrichment**: *true, false*
 - **do_saturation_curve**: *true, false*
 - **do_raw_peak_annotation**: *true, false*
 - **do_diffbind_peak_annotatio**: *true, false*
 - **do_gene_set_enrichment**: *true, false*
 - **do_bigwig**: *true, false*
 - **do_chromatin_state**: *true, false*

## Experiment
 - **specie**: *worm, fly, mouse or human*
 - **experiment_types**: *both, atac or mRNA*.  
To indicate if the run should analyze ATAC data only, mRNA data only, or both kind of data.
- **use_input_control**: *false, true*.  
Should a gDNA input control be used for ATAC-Seq analysis to remove [greylist regions](https://rdrr.io/bioc/DiffBind/man/dba.blacklist.html) with DiffBind.

## ATAC reads analysis
- **sam_MAPQ_threshold**: *30*. MAPQ threshold for filtering low quality reads.
- **binsize_bigwig_creation**: *10*. Size of the bins for the creation of the bigwig files with DeepTools. Larger numbers will make the analysis faster but the resolution smaller. Bins of 10000 can be used for a fast analysis.
- **binsize_bigwig_correlation**: *10000*. Size of the bins for the correlation analysis with DeepTools.
- **nb_1bp_site_to_sample_for_coverage**: *10000*. Number of sites to sample for the coverage plots with DeepTools.
- **nb_sampled_reads**: *1000000*. Number of reads to sample for the resampling analysis.
- **macs2_qvalue**: *5e-2*. Number of reads to sample for the resampling analysis.
- **macs2_mappable_genome_size**: *9e7*. Number of reads to sample for the resampling analysis.


//// genes quantification
fragment_len  = '180'
fragment_sd   = '20'
bootstrap     = '100'

//// promoters annotation
promoter_up_macs2_peaks      = 1500
promoter_down_macs2_peaks    = 500
promoter_up_diffbind_peaks   = 1500
promoter_down_diffbind_peaks = 500

//// volcano and MA plots
fdr_threshold_diffbind_plots = 0.05
fdr_threshold_sleuth_plots   = 0.05

//// splitting DA results in subset
threshold_type_for_splitting_subsets   = 'FDR'  // options: 'FDR' or 'rank'
threshold_values_for_splitting_subsets = [ 1.3 ] // rank cuttoff or -log10 FDR thresholds (here: -log10(0.05) ~= 1.3)
fold_changes_for_splitting_subsets     = [ 'up', 'down' ] // For now this argument should not be touched. But an 'all' argument may appear in the future
peak_assignment_for_splitting_subsets  = [ 'all', 'prom', 'distNC' ] // this argument can be any of: 'all', '8kb', '3kb', '2u1d', 'TSS', 'genProm', 'genic', 'prom' or 'distNC'

//// filtering DA results with too few entries
min_entries_DA_bed        = 2
min_entries_DA_genes_sets = 1

//// functional annotation databases; possible values: BP, CC, MF, KEGG
params.func_anno_databases = ['BP', 'KEGG']
params.use_nda_as_bg_for_func_anno = 'F'

// motifs test can be eitheir hypergeometric (default, any value) or binomial
motifs_test_type = 'binomial'

// threshold for saving enrichment results (so 1 allows to save all results)
fdr_filter_tables =       [  1,    1, 1,     1,             1,            1,          1,              1,    1,             1 ]
fdr_filter_tables_names = [ 'mRNA_detailed', 'ATAC_detailed', 'res_simple', 'res_filter', 'func_anno', 'genes_self', 'peaks_self', 'chrom_states', 'CHIP', 'motifs' ]

//// heatmaps and barplots
add_var_to_plot         = 'L2OR' // can be either: 'Log2OddRatios', 'loglog_pvalue' or 'DA_overlap_count'
threshold_plot_adj_pval = 0.05   // at least one enrichment should have a pvalue lower than that otherwise no plot will be made
up_down_pattern         = 'UUDD' //  can be either: UDUD (up, down, up, down...) or UUDD (up, up, ..., down, down ...) 

// table output: csv or excel
// params.tables_output_type = 'excel'
params.save_tables_as_csv   = false
params.save_tables_as_excel = true
// params.save_merged_tables = true
// params.save_individual_tables = true // => these two parameters should be implemented later

params.excel_add_conditional_formatting = 'T'
params.excel_max_width = 40


nb_threads_pigz       = 6
nb_threads_botwie2    = 6
nb_threads_deeptools  = 6
nb_threads_kallisto   = 6
nb_threads_homer      = 6