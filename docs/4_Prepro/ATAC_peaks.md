

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Plots](/docs/6_Enrich/Plots.md), [Reports](/docs/6_Enrich/Reports.md)

[](END_OF_MENU)


# List of processes

- [Preprocessing](#Preprocessing)
  - [ATAC_peaks__calling_peaks](#ATAC_peaks__calling_peaks)
  - [ATAC_peaks__splitting_multi_summits_peaks](#ATAC_peaks__splitting_multi_summits_peaks)
  - [ATAC_peaks__removing_blacklisted_regions](#ATAC_peaks__removing_blacklisted_regions)
  - [ATAC_peaks__removing_input_control_peaks](#ATAC_peaks__removing_input_control_peaks)
  - [ATAC_peaks__removing_specific_regions](#ATAC_peaks__removing_specific_regions)

- [Quality Controls](#Quality-Controls)
  - [ATAC_QC_peaks__computing_and_plotting_saturation_curve](#ATAC_QC_peaks__computing_and_plotting_saturation_curve)
  - [ATAC_QC_peaks__annotating_macs2_peaks](#ATAC_QC_peaks__annotating_macs2_peaks)
  - [ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_each_sample](#ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_each_sample)
  - [ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_all_samples_grouped](#ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_all_samples_grouped)


# Preprocessing

## ATAC_peaks__calling_peaks

### Description
Inputs are reads in bed files that are 1 base pair long (the 5' end), and that have been adjusted for the [shift of the transposase](https://doi.org/10.1038/nmeth.2688).  
Peaks are called with [MACS2](https://doi.org/10.1186/gb-2008-9-9-r137). The `macs2 callpeak` function is called with these arguments: `-f BED -nomodel --shift -75 --extsize 150 --call-summits`.  
The `--call-summits` argument allows to call multiple summits for each peaks. This allows to split peaks in the next process.  
The `--extsize 150` argument allows to extend the reads by 150 bp in the 3' direction. This way fragments are centered in the original 5' location.  
The combination of the `--shift` and `--extsize` arguments are often recommended for ATAC-Seq data (see these links for references and explanations: [1](https://groups.google.com/g/macs-announcement/c/4OCE59gkpKY/m/v9Tnh9jWriUJ), [2](https://github.com/macs3-project/MACS/issues/145#issuecomment-742593158), [3](https://github.com/macs3-project/MACS/discussions/435), [4](https://twitter.com/XiChenUoM/status/1336658454866325506)).  
We use a size of 150 base pair as it is approximately the size of a nucleosome. 

>**_Note_:** the two ends of each read pair are analyzed separately as they both provide valuable and separate information. See [here](https://twitter.com/XiChenUoM/status/1336658454866325506).


### Parameters
- **_params.macs2_qvalue_**: q-value (minimum FDR) cutoff to call significant regions.

### Outputs
- **Raw peaks** (`${id}__macs2_peaks.narrowPeak`) if **_params.save_bed_type = 'all'_** in `Processed_Data/1_Preprocessing/ATAC__peaks__raw`.


## ATAC_peaks__splitting_multi_summits_peaks

### Description
MACS2 peaks with multiple summits are split, with a boundary set in the middle of neighboring summits.  
The script from this process was written by [Aaron C Daugherty](https://github.com/brunetlab/CelegansATACseq/blob/master/Fig1/splitMACS2SubPeaks.pl) for the [first ATAC-Seq paper in *C. elegans*](http://www.genome.org/cgi/doi/10.1101/gr.226233.117).

### Outputs
- **Split peaks** (`${id}__split_peaks.narrowPeak`) if **_params.save_bed_type = 'all'_** in `Processed_Data/1_Preprocessing/ATAC__peaks__split`.


## ATAC_peaks__removing_blacklisted_regions

### Description
Any peak that has any overlap with a blacklisted region is discarded.

### Outputs
- **Kept and discarded peaks** (`${id}__peaks_kept_after_blacklist_removal.bed` and `${id}__peaks_lost_after_blacklist_removal.bed`) if **_params.save_bed_type = 'all'_** in `Processed_Data/1_Preprocessing/ATAC__peaks__split__no_BL`.



## ATAC_peaks__removing_input_control_peaks

### Description
If an input control is included in the experiment, and **_params.use_input_control_** is true, then peaks overlapping with input control peaks are removed.

### Parameterss
- **_params.input_control_overlap_portion_**: threshold of the fraction of overlapping input control peaks to remove peaks. The percentage is regarding the treatment/sample peaks, not the input control peaks. Default: 0.2.

### Outputs
- **Kept and discarded peaks** (`${id}__peaks_kept_after_input_control_removal.bed` and `${id}__peaks_lost_after_input_control_removal.bed`) if **_params.save_bed_type = 'all'_** in `Processed_Data/1_Preprocessing/ATAC__peaks__split__no_BL_input`.


## ATAC_peaks__removing_specific_regions

### Description
This process takes as input peaks from a comparison and remove any peaks that are in a region to remove for any of the two sample_id to compare. This is mostly useful for RNAi experiments that induce a large signal in ATAC-Seq, but it can be used to remove any arbitrary region.  

>**_Note_:** the reason why this process is here and not upstream (before aggregating replicates and comparisons) is because we want to remove in all bed files the peaks that are in specific regions (i.e. RNAi) that we want to avoid. This is because, Diffbind (i.e. the Differential Binding Analysis (DBA) tool) will consider all peaks for his analysis (i.e. differential abundance for ATAC-Seq), so if we remove one such peak in just one of the two samples to compare, if it is still present in the other sample then it will be included in the analysis and it will likely be found as differential bound during the DBA. e.g.: let's say we compare daf-16 RNAi vs control. If there is a MACS2 peak at daf-16 in one of the control condition's replicate, then even if we remove this peak in the daf-16 RNAi condition's replicates, it will still be included in the final analysis.

### Parameterss
- **_params.design__regions_to_remove**: path to the file containing the regions to remove (see the [Design](/docs/3_Inputs/Design.md) section for details). Default: 'Design/regions_to_remove.tsv'.

### Outputs
- **Kept and discarded peaks** (`${id}__peaks_kept_after_specific_regions_removal.bed` and `${id}__peaks_lost_after_specific_regions_removal.bed`) if **_params.save_bed_type = 'all'_** in `Processed_Data/1_Preprocessing/ATAC__peaks__split__no_BL_input`.






# Quality Controls

## ATAC_QC_peaks__computing_and_plotting_saturation_curve

### Description
A saturation curve is made to help assessing if samples have been sequenced deep enough (see [here](https://doi.org/10.1038/nrg3642) for details on saturation curve analysis).  
MACS2 is used to sample reads (1 base pair reads in bed format) from 10% to 100%, by 10% increase.  
Then, peaks are called with MACS2 for each sample. 
Finally, a plot is made in R showing the number of peaks (y-axis) by sequencing depth (x-axis). 

### Parameterss
- **_params.do_saturation_curve_**: enable or disable this process. Default: true.

### Outputs
- **Saturation curves** (`${id}__saturation_curve.pdf`) in `Figures_Individual/1_Preprocessing/ATAC__peaks__saturation_curve`


## ATAC_QC_peaks__annotating_macs2_peaks

### Description
Peaks are annotated with [ChIPseeker](http://dx.doi.org/10.1093/bioinformatics/btv145). Each peak is assigned to its closest gene using the [annotatePeak function](https://github.com/YuLab-SMU/ChIPseeker/blob/master/R/annotatePeak.R).

### Parameters
- **_params.do_raw_peak_annotation_**: to enable or disable this process. Default: true.
- **_params.promoter_up_macs2_peaks_**: promoter start; upstream from TSS site.
- **_params.promoter_down_macs2_peaks_**: promoter end; downstream from TSS site.

### Outputs
- **Annotated peaks R objects** (`${id}__annotated_peaks.pdf`) in `Processed_Data/1_Preprocessing/ATAC__peaks__annotated_rds`.


## ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_each_sample

### Description
Using ChIPseeker and [ggplot2](https://ggplot2.tidyverse.org/) to plot coverage and average profile around TSS for each sample (one plot type per sample).

### Parameters
- **_params.promoter_up_macs2_peaks_**: promoter start; upstream from TSS site.
- **_params.promoter_down_macs2_peaks_**: promoter end; downstream from TSS site.

### Outputs
- **Coverage plots** (`${id}__coverage.pdf`) in `Figures_Individual/1_Preprocessing/ATAC__peaks__coverage`.
- **Average profile plots** (`${id}__average_profile.pdf`) in `Figures_Individual/1_Preprocessing/ATAC__peaks__average_profile`.
  

## ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_all_samples_grouped

### Description
Using ChIPseeker and ggplot2 to plot coverage and average profile around TSS for all samples grouped (one plot type for all samples).

### Parameters
- **_params.promoter_up_macs2_peaks_**: promoter start; upstream from TSS site.
- **_params.promoter_down_macs2_peaks_**: promoter end; downstream from TSS site.

### Outputs
- **Average profile plots** (`ATAC__peaks__average_profile.pdf`) made by the [plotAvgProf](https://rdrr.io/bioc/ChIPseeker/man/plotAvgProf.html) function
- **Annotation barplots** (`ATAC__peaks__annotation_barplot.pdf`) made by the [plotAnnoBar](https://rdrr.io/bioc/ChIPseeker/man/plotAnnoBar.data.frame.html) function
- **Distance to TSS** (`ATAC__peaks__distance_to_TSS.pdf`)  made by the [plotDistToTSS](https://rdrr.io/bioc/ChIPseeker/man/plotDistToTSS.data.frame.html) function
  - in `Figures_Individual/1_Preprocessing/ATAC__peaks__grouped_plots`.


