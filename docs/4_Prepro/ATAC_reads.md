

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [DBA](/docs/5_DA/DBA.md), [DGEA](/docs/5_DA/DGEA.md), [Split](/docs/5_DA/Split.md), [Outputs](/docs/5_DA/Outputs.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Outputs](/docs/6_Enrich/Outputs.md)

[](END_OF_MENU)


## Processes

 - ATAC_reads__merging_reads: samples sequenced multiple times (that have the same id) are merged
 - ATAC_reads__trimming_reads: adapters are removed skewer 
 - ATAC_reads__aligning_reads
 - ATAC_reads__removing_low_quality_reads
 - ATAC_reads__marking_duplicated_reads
 - ATAC_reads__removing_duplicated_reads
 - ATAC_reads__removing_reads_in_mitochondria_and_small_contigs
 - ATAC_reads__converting_bam_to_bed_and_adjusting_for_Tn5
 
 ATAC_QC_reads__running_fastqc
 ATAC_QC_reads__computing_bigwig_tracks_and_plotting_coverage
 ATAC_QC_reads__computing_and_plotting_bigwig_tracks_correlations
 ATAC_QC_reads__plotting_insert_size_distribution
 ATAC_QC_reads__sampling_aligned_reads
 ATAC_QC_reads__measuring_overlap_with_genomic_regions
 ATAC_QC_reads__estimating_library_complexity
 ATAC_QC_reads__sampling_trimmed_reads
 ATAC_QC_reads__aligning_sampled_reads
 ATAC_QC_reads__gathering_all_stat
 ATAC_QC_reads__gathering_all_samples
 ATAC_QC_reads__splitting_stat_for_multiqc
 ATAC_QC_reads__running_multiQC


## steps

If there are duplicate 


## parameters

## outputs






