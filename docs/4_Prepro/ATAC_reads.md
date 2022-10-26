

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


# List of processes

- [Preprocessing](#Preprocessing)
  - [ATAC_reads__merging_reads](#ATAC_reads__merging_reads)
  - [ATAC_reads__trimming_reads](#ATAC_reads__trimming_reads)
  - [ATAC_reads__aligning_reads](#ATAC_reads__aligning_reads)
  - [ATAC_reads__removing_low_quality_reads](#ATAC_reads__removing_low_quality_reads)
  - [ATAC_reads__marking_duplicated_reads](#ATAC_reads__marking_duplicated_reads)
  - [ATAC_reads__removing_duplicated_reads](#ATAC_reads__removing_duplicated_reads)
  - [ATAC_reads__removing_reads_in_mitochondria_and_small_contigs](#ATAC_reads__removing_reads_in_mitochondria_and_small_contigs)
  - [ATAC_reads__converting_bam_to_bed_and_adjusting_for_Tn5](#ATAC_reads__converting_bam_to_bed_and_adjusting_for_Tn5)  


- [Quality Controls](#Quality-Controls)
  - [ATAC_QC_reads__computing_bigwig_tracks_and_plotting_coverage](#ATAC_QC_reads__computing_bigwig_tracks_and_plotting_coverage)
  - [ATAC_QC_reads__computing_and_plotting_bigwig_tracks_correlations](#ATAC_QC_reads__computing_and_plotting_bigwig_tracks_correlations)
  - [ATAC_QC_reads__plotting_insert_size_distribution](#ATAC_QC_reads__plotting_insert_size_distribution)
  - [ATAC_QC_reads__sampling_aligned_reads](#ATAC_QC_reads__sampling_aligned_reads)
  - [ATAC_QC_reads__measuring_overlap_with_genomic_regions](#ATAC_QC_reads__measuring_overlap_with_genomic_regions)
  - [ATAC_QC_reads__estimating_library_complexity](#ATAC_QC_reads__estimating_library_complexity)
  - [ATAC_QC_reads__sampling_trimmed_reads](#ATAC_QC_reads__sampling_trimmed_reads)
  - [ATAC_QC_reads__aligning_sampled_reads](#ATAC_QC_reads__aligning_sampled_reads)
  - [ATAC_QC_reads__gathering_all_stat](#ATAC_QC_reads__gathering_all_stat)
  - [ATAC_QC_reads__gathering_all_samples](#ATAC_QC_reads__gathering_all_samples)
  - [ATAC_QC_reads__splitting_stat_for_multiqc](#ATAC_QC_reads__splitting_stat_for_multiqc)


# Preprocessing

## ATAC_reads__merging_reads

### Description
Samples sequenced multiple times (that have the same id) are merged.

### Outputs
- **Merged reads**: `Processed_Data/1_Preprocessing/ATAC__reads__fastq_merged/${sample}_R{1,2}_merged.fastq.gz` if *params.save_fastq_type = 'all'* in .


## ATAC_reads__trimming_reads
 
### Description
ATAC-Seq adaptors are trimmed using [Skewer](https://doi.org/10.1186/1471-2105-15-182), then they are compressed in parallel by [PIGZ](https://zlib.net/pigz/).

### Parameters
- **_params.pigz__nb_threads_**: number of threads used for parallel compression. Default: 6.

### Outputs
- **Log files**: 
  - `${sample}_skewer_trimming.log`
  - `${sample}_pigz_compression.log`
- **Trimmed reads**: `*_R{1,2}_trimmed.fastq` if *params.save_fastq_type = 'all'*

### Output folders
- `Processed_Data/1_Preprocessing/ATAC__reads__fastq_trimmed`.


## ATAC_reads__aligning_reads

### Description
Reads are aligned to the reference genome by [Bowtie2](https://doi.org/10.1038/nmeth.1923) and the resulting SAM files are converted to BAM with [SAMtools](https://pubmed.ncbi.nlm.nih.gov/19505943/). SAMtools is also used to count the number of aligned reads per category.

### Parameters
- **_params.bowtie2__nb_threads_**: number of threads used by Bowtie2. Default: 6.

### Outputs
- **Bowtie 2 alignment metrics**: `${sample}_bowtie2_align_metrics.txt`
- **Number of aligned reads per category**: `${sample}.qc`
- **Aligned reads**: `${sample}.bam` if *params.save_bam_type = 'all'* 

### Output folders
- `Processed_Data/1_Preprocessing/ATAC__reads__fastq_trimmed`.


## ATAC_reads__removing_low_quality_reads

### Description
Low quality reads with these attributes are filtered: unmapped, mate unmapped, no primary alignment, low MAPQ (quality score). Reads are sorted and the number of aligned reads per category is determined with SAMtools.

### Parameters
- **_params.sam_MAPQ_threshold_**: MAPQ threshold. Default: 30.

### Outputs
- **Number of aligned reads per category**: `${sample}__filter_LQ.qc`
- **Aligned reads**: `${sample}__filter_LQ.bam` if *params.save_bam_type = 'all'* 

### Output folders
- `Processed_Data/1_Preprocessing/ATAC__reads__bam_no_lowQ`.


## ATAC_reads__marking_duplicated_reads

### Description
[Picard](https://broadinstitute.github.io/picard/) is used to mark duplicated reads.

### Parameters
- **_params.memory_picard_**: maximum memory used by Picard. Default: '20G'.


## ATAC_reads__removing_duplicated_reads

### Description
Duplicated reads are removed with SAMtools, an index is build, and the number of aligned reads per category is determined with SAMtools.

### Outputs
- **Number of aligned reads per category**: `${sample}__dup_rem.qc`
- **Aligned reads**: `${sample}__dup_rem.bam` if *params.save_bam_type = 'all'* 

### Output folders
- `Processed_Data/1_Preprocessing/ATAC__reads__bam_no_lowQ_dupli`.


## ATAC_reads__removing_reads_in_mitochondria_and_small_contigs

### Description
Samtools is used to remove reads mapping to the mitochondrial chromosome or to small contigs. This is defined as regions (contigs and chromosomes) with less than 5 genes (see [References](/docs/2_Install/References.md#Pipeline-to-get-references)). Number of aligned reads per chromosomes and per category are also computed with SAMtools.

### Outputs
- **Number of aligned reads per category**: `${sample}__no_mito.qc`
- **Number of reads per chromosome before and after filtering**: .txt file)
- **Aligned reads**: `${sample}__no_mito.bam` if *params.save_bam_type = 'all' or 'last'* 

### Output folders
- `Processed_Data/1_Preprocessing/ATAC__reads__bam_no_lowQ_dupli_mito`.


 
## ATAC_reads__converting_bam_to_bed_and_adjusting_for_Tn5

### Description
This process converts bam files to bed files, and use the insert length as the score. It then adjusts for the shift of the transposase to account for the [9 bp offset of the transposase](https://doi.org/10.1038/nmeth.2688). This is done by shifting reads on the + strand by +4 base pairs and reads on the - strand by -5 base pairs. Finally, it keeps only the 5' end of each reads (and thus each read becomes 1 base pair long), and create a sorted and indexed bam file from the final adjusted bed file. The bam file is sent to [DiffBind](https://doi.org/10.1038/nature10730)] for Differential Binding analysis. The bed file is sent for custom quality controls processes (see below), for computing and plotting saturation curve, and for calling macs2 peaks.

### Outputs
  - **Aligned reads**: `Processed_Data/1_Preprocessing/ATAC__reads__bam_asBed_atacShift/${sample}__1bp_shifted_reads.bam` if *params.save_1bp_bam = true*.


 
# Quality Controls


## ATAC_QC_reads__running_fastqc

### Description
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is runned to perform standard quality controls checks on sequenced data.

### Parameters
- **_params.fastqc__nb_threads_**: number of threads used by FastQC. Default: 2.

### Outputs
- **Reads quality control reports**: `*_R{1,2}_fastqc.html`

## Output folders
 - `Processed_Data/1_Preprocessing/ATAC__reads__fastqc_raw` for raw reads
 - `Processed_Data/1_Preprocessing/ATAC__reads__fastqc_trimmed` for trimmed reads.


## ATAC_QC_reads__computing_bigwig_tracks_and_plotting_coverage

### Description
Coverage profiles (i.e. bigwig files) and plots are generated by [DeepTools](https://doi.org/10.1093/nar/gkw257).

### Parameters
- **_params.do_bigwig_**: enable or disable this process. Default: true.
- **_params.deeptools__binsize_bigwig_creation_**: size of the bins in the bigwig file. Smaller values increase computation time. Default: 10000.
- **_params.deeptools__nb_threads_**: number of threads used by DeepTools. Default: 6.
- **_params.deeptools__nb_of_1_bp_samples_**: number of 1 bp sites to sample for the coverage plots. Default: 10000.
- **_params.deeptools__normalization_method_**: normalization method to use when creating BigWig files. See [here](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html) for options. Default: 'none'.

### Outputs
- **[BigWig files](https://deeptools.readthedocs.io/en/latest/content/tools/bamCoverage.html)**: `Processed_Data/1_Preprocessing/ATAC__reads__bigwig_raw/${sample}.bw`.
- **[Coverage plots](https://deeptools.readthedocs.io/en/latest/content/tools/plotCoverage.html)**: 
  - `Figures_Individual/1_Preprocessing/ATAC__reads__coverage/${sample}__reads_coverage.pdf`
  - `Figures_Merged/1_Preprocessing/ATAC__reads__coverage.pdf`
![Reads Coverage](/docs/examples/png/ctl_1__reads_coverage.png)


## ATAC_QC_reads__computing_and_plotting_bigwig_tracks_correlations

### Description
First, the coverage matrix of all samples is computed with [DeepTools](https://doi.org/10.1093/nar/gkw257)'s function multiBigwigSummary. Then PCA and Correlations plots are generated also with DeepTools, while excluding blacklisted regions.
The PCA is computed with filtering genes to keep only the top (100, 1000, 5000) most variables rows in the matrix, and with or without genomic DNA control (i.e. input) if present.  
Correlation matrix is computed for spearman and pearson correlation, with or without outliers, and with or without genomic DNA control (i.e. input) if present.
If an input control is included (genomic DNA), then two sets of plots will be made: with or without the input.

### Parameters
- **_params.deeptools__binsize_bigwig_creation_**: size of the bins in the coverage matrix. Smaller values increase computation time. Default: 10000.

### Outputs
- **[PCA plots](https://deeptools.readthedocs.io/en/latest/content/tools/plotPCA.html)**: 
  - `Figures_Individual/1_Preprocessing/ATAC__reads__PCA/pca_top{100,1000,5000}_${gDNA_PRESENT}_pca.pdf`
  - `Figures_Merged/1_Preprocessing/ATAC__reads__PCA.pdf`.
![Reads PCA](/docs/examples/png/pca_top5000_without_control_pca.png)
- **[Correlation plots](https://deeptools.readthedocs.io/en/latest/content/tools/plotCorrelation.html)**: 
  - `Figures_Individual/1_Preprocessing/ATAC__reads__correlations/{pearson,spearman}_correlation_heatmap_{with,without}_outliers_${gDNA_PRESENT}_cor.pdf`
  - `Figures_Merged/1_Preprocessing/ATAC__reads__correlations.pdf`.
![Reads PCA](/docs/examples/png/spearman_correlation_heatmap_without_outliers_without_control_cor.png)

  
## ATAC_QC_reads__plotting_insert_size_distribution

### Description
Insert size plots are made with [Picard](https://broadinstitute.github.io/picard/).

### Parameters
- **_params.memory_picard_**: maximum memory used by Picard. Default: '20G'.

### Outputs
- **[insert size plots](https://gatk.broadinstitute.org/hc/en-us/articles/360037225252-CollectInsertSizeMetrics-Picard-)**: 
  - `Figures_Individual/1_Preprocessing/ATAC__reads__insert_size/${sample}__insert_size.pdf` 
  - `Figures_Merged/1_Preprocessing/ATAC__reads__insert_size.pdf`.
![Insert size](/docs/examples/png/ctl_1__insert_size.png)

 
## ATAC_QC_reads__sampling_aligned_reads

### Description
This process takes as input non-filtered aligned reads (coming straight out of Bowtie2).  
These are filterd by keeping only mapped reads (flag 0x4) that are not secondary (flag 0x100, multimappers) or supplementary (flag 0x800, chimeric entries) and transformed to the sam format with SAMtools. This sam file is sent to [estimating_library_complexity](/docs/4_Prepro/ATAC_reads.md#ATAC_QC_reads__estimating_library_complexity).  
Then a certain number of reads are randomly sampled and a bam file is created, sorted and indexed. This bam file is sent to [measuring_overlap_with_genomic_regions](/docs/4_Prepro/ATAC_reads.md#ATAC_QC_reads__measuring_overlap_with_genomic_regions) and [gathering_all_stat](/docs/4_Prepro/ATAC_reads.md#ATAC_QC_reads__gathering_all_stat).  
Finally, SAMtools is used again to determine the number of total pairs and the number of aligned pairs in the original, non-sampled bam file. This information is sent to [gathering_all_stat](/docs/4_Prepro/ATAC_reads.md#ATAC_QC_reads__gathering_all_stat).

### Parameters
- **_params.nb_sampled_aligned_reads_**: Number of aligned reads to sample. Default: 1000000.
 
 
## ATAC_QC_reads__measuring_overlap_with_genomic_regions

### Description
The overlap of sampled aligned reads with various annotated regions (exons, genes, intergenic regions, introns and promoters) is computed with [BEDTools](https://doi.org/10.1093/bioinformatics/btq033).

 
## ATAC_QC_reads__estimating_library_complexity

### Description
[Picard](https://broadinstitute.github.io/picard/) is used to estimate library complexity of our set of sampled aligned reads.

### Parameters
- **_params.memory_picard_**: maximum memory used by Picard. Default: '20G'.


## ATAC_QC_reads__sampling_trimmed_reads

### Description
[BBMap](https://sourceforge.net/projects/bbmap/) is used to sample trimmed reads. 

### Parameters
- **_params.nb_sampled_trimmed_reads_**: Number of trimmed reads to sample. Default: 1000000.


## ATAC_QC_reads__aligning_sampled_reads

### Description
Sampled trimmed reads are aligned to the reference genome under study (worm, fly, mouse or human) and to the contaminant specie's reference genome (i.e. *E. coli* OP50 strain for all species currently). Then, the number of aligned reads per category is determined with SAMtools for both organisms.

### Parameters
- **_params.botwie2__nb_threads_**: number of threads used by Bowtie2. Default: 6.

 
## ATAC_QC_reads__gathering_all_stat

### Description
This process gather all reads statistics from previous processes (starting with *ATAC_QC_reads__*) for a given sample, computes a few more stastistics, and creates a summary file.  

Here is a description/summary of how the reported statistics were generated:  

- Number of paired end reads:
  - Raw Pairs: number of total reads in the original bam file (after Bowtie2 alignment, without any filtering step)
  - Aligned Pairs: number of aligned reads in the original bam file (after Bowtie2 alignment, without any filtering step)
  - Final Pairs: number of reads in the final bed file (1 base pair reads after all filtering steps)

- Percentage of sampled trimmed reads that aligned to: 
  - the reference genome under study
  - the contaminant specie's genome.

- Sampled aligned reads are used to estimate:
  - The percentage of read that overlap with the mitochondrial genome and with genomic annotations (promoters, exons, introns, intergenic regions and genic regions)
  - Library complexity


## ATAC_QC_reads__gathering_all_samples

### Description
Reads statistics of all samples are combined into a single file.

 
## ATAC_QC_reads__splitting_stat_for_multiqc

### Description
Variable names are formated in R, and preprocessing is made for making a file suitable for [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/).

 
## ATAC_QC_reads__running_multiQC

### Description
A MultiQC html report is made that aggregates all basic FastQC quality control files, and the custom statistics files generated by cactus.

### Outputs
- **MultiQC report**: `ATAC__multiQC.html`.
 
### Output folders
- `Figures_Individual/1_Preprocessing`
- `Figures_Merged/1_Preprocessing`
