

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Plots](/docs/6_Enrich/Plots.md), [Reports](/docs/6_Enrich/Reports.md)

[](END_OF_MENU)



## ATAC_reads__merging_reads

  ### Description
Samples sequenced multiple times (that have the same id) are merged.
  
  ### Outputs
#### Files
 - **Fasqt files after merging** if *params.save_fastq_type = 'all'*
 #### Folder
*Processed_Data/1_Preprocessing/ATAC__reads__fastq_merged*


## ATAC_reads__trimming_reads
 
  ### Description
ATAC-Seq adaptors are trimmed using (Skewer)[https://doi.org/10.1186/1471-2105-15-182], then they are compressed in parallel by 
[PIGZ](https://zlib.net/pigz/).

  ### Parameters
  - *params.nb_threads_pigz* controls the number of threads used for parallel compression.

  ### Outputs
  #### Files
  - **Fasqt files after trimming** if *params.save_fastq_type = 'all'*
  - **Trimming and compression log files**
  #### Folder
  *Processed_Data/1_Preprocessing/ATAC__reads__fastq_trimmed*

## ATAC_reads__aligning_reads

 ### Description

 ### Parameters

 ### Outputs
 #### Files
 - **Bowt**: (${id}_skewer_trimming.log and ${id}_pigz_compression.log) 
 - **Fasqt files after trimming**: (${id}_R1_trim.fastq.gz and ${id}_R2_trim.fastq.gz) if *params.save_bam_type = 'all'* 
 #### Folder
 *Processed_Data/1_Preprocessing/ATAC__reads__fastq_trimmed*

## ATAC_reads__removing_low_quality_reads

 ### Description

 ### Parameters

 ### Outputs


## ATAC_reads__marking_duplicated_reads

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_reads__removing_duplicated_reads

 ### Description

 ### Parameters

 ### Outputs
 
## ATAC_reads__removing_reads_in_mitochondria_and_small_contigs

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_reads__converting_bam_to_bed_and_adjusting_for_Tn5

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__running_fastqc

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__computing_bigwig_tracks_and_plotting_coverage

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__computing_and_plotting_bigwig_tracks_correlations

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__plotting_insert_size_distribution

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__sampling_aligned_reads

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__measuring_overlap_with_genomic_regions

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__estimating_library_complexity

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__sampling_trimmed_reads

 ### Description

 ### Parameters

 ### Outputs
 

## ATAC_QC_reads__aligning_sampled_reads

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__gathering_all_stat

 ### Description

 ### Parameters

 ### Outputs
 

## ATAC_QC_reads__gathering_all_samples

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__splitting_stat_for_multiqc

 ### Description

 ### Parameters

 ### Outputs
 
 
## ATAC_QC_reads__running_multiQC

 ### Description

 ### Parameters

 ### Outputs
 


