

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Plots](/docs/6_Enrich/Plots.md), [Reports](/docs/6_Enrich/Reports.md)

[](END_OF_MENU)

[Analysis](#Analysis)
 - [ATAC_peaks__calling_peaks](#ATAC_peaks__calling_peaks)
 - [ATAC_peaks__splitting_multi_summits_peaks](#ATAC_peaks__splitting_multi_summits_peaks)
 - [ATAC_peaks__removing_blacklisted_regions](#ATAC_peaks__removing_blacklisted_regions)
 - [ATAC_peaks__removing_input_control_peaks](#ATAC_peaks__removing_input_control_peaks)
 - [ATAC_peaks__removing_specific_regions](#ATAC_peaks__removing_specific_regions)

[Quality Controls](#Quality-Controls)
 - [ATAC_QC_peaks__computing_and_plotting_saturation_curve](#ATAC_QC_peaks__computing_and_plotting_saturation_curve)
 - [ATAC_QC_peaks__annotating_macs2_peaks](#ATAC_QC_peaks__annotating_macs2_peaks)
 - [ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_each_sample](#ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_each_sample)
 - [ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_all_samples_grouped](#ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_all_samples_grouped)


# Analysis

## ATAC_peaks__calling_peaks

### Description
Inputs are reads in bed files that are 1 base pair long (the 5' end), and that have been adjusted for the [shift of the transposase](https://doi.org/10.1038/nmeth.2688).
Peaks are called with [MACS2](https://doi.org/10.1186/gb-2008-9-9-r137). These options are specified: `-f BED -nomodel --shift -75 --extsize 150`.  
The argument `--extsize 150` indicates to extend the reads by 150 bp in the 3' direction. This way fragments are centered in the original 5' location. 
The combination of the --shift and --extsize arguments are often recommended for ATAC-Seq data (see these links for references and explanations: [1](https://groups.google.com/g/macs-announcement/c/4OCE59gkpKY/m/v9Tnh9jWriUJ), [2](https://github.com/macs3-project/MACS/issues/145#issuecomment-742593158), [3](https://github.com/macs3-project/MACS/discussions/435), [4](https://twitter.com/XiChenUoM/status/1336658454866325506)).
We use a size of 150 base pair as it is approximately the size of a nucleosome. 

> Note: the two ends of each read pair are analyzed separately as they both provide valuable and separate information. See [here](https://twitter.com/XiChenUoM/status/1336658454866325506).


### Parameters
- **_params.SDS**: EER.

### Outputs
- **SS** (.log files)
- **FF** if **_params.XX = 'EE'_**
  - in `XX`.


## ATAC_peaks__splitting_multi_summits_peaks

### Description

### Parameters
- **_params.SDS**: EER.

### Outputs
- **SS** (.log files)
- **FF** if **_params.XX = 'EE'_**
  - in `XX`.


## ATAC_peaks__removing_blacklisted_regions

### Description

### Parameters
- **_params.SDS**: EER.

### Outputs
- **SS** (.log files)
- **FF** if **_params.XX = 'EE'_**
  - in `XX`.


## ATAC_peaks__removing_input_control_peaks

### Description

### Parameters
- **_params.SDS**: EER.

### Outputs
- **SS** (.log files)
- **FF** if **_params.XX = 'EE'_**
  - in `XX`.


## ATAC_peaks__removing_specific_regions

### Description

### Parameters
- **_params.SDS**: EER.

### Outputs
- **SS** (.log files)
- **FF** if **_params.XX = 'EE'_**
  - in `XX`.





# Quality Controls

## ATAC_QC_peaks__computing_and_plotting_saturation_curve

### Description

### Parameters
- **_params.SDS**: EER.

### Outputs
- **SS** (.log files)
- **FF** if **_params.XX = 'EE'_**
  - in `XX`.


## ATAC_QC_peaks__annotating_macs2_peaks

### Description

### Parameters
- **_params.SDS**: EER.

### Outputs
- **SS** (.log files)
- **FF** if **_params.XX = 'EE'_**
  - in `XX`.


## ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_each_sample

### Description

### Parameters
- **_params.SDS**: EER.

### Outputs
- **SS** (.log files)
- **FF** if **_params.XX = 'EE'_**
  - in `XX`.


## ATAC_QC_peaks__plotting_annotated_macs2_peaks_for_all_samples_grouped

### Description

### Parameters
- **_params.SDS**: EER.

### Outputs
- **SS** (.log files)
- **FF** if **_params.XX = 'EE'_**
  - in `XX`.


