

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Plots](/docs/6_Enrich/Plots.md), [Reports](/docs/6_Enrich/Reports.md)

[](END_OF_MENU)


# List of processes

- [Preprocessing](#Preprocessing)
  - [MRNA__quantifying_transcripts_abundances](#MRNA__quantifying_transcripts_abundances)

- [Quality Controls](#Quality-Controls)
  - [MRNA_QC__running_fastqc](#MRNA_QC__running_fastqc)
  - [MRNA_QC__running_MultiQC](#MRNA_QC__running_MultiQC)



# Preprocessing

 
## MRNA__quantifying_transcripts_abundances

### Description
Quantification of transcripts abundance is made with [kallisto](https://doi.org/10.1038/nbt.3519) in single or paired-end mode according to the input type. Briefly, Kallisto is based on the recent idea of a pseudoalignment that allows fasta and precise alignment-free quantification of transcripts (more details [here](https://pachterlab.github.io/kallisto/about)).

### Parameterss
- **_params.nb_threads_kallisto_**: number of threads used by kallisto. Default: 6.
- **_params.bootstrap_**: Number of bootstrap samples. Default: '100'.
- **_params.fragment_len_**: Estimated average fragment length. For single end only. Default: '180'.
- **_params.fragment_sd_**: Estimated standard deviation of fragment length. For single end only. Default: '20'.
- **_params.nb_threads_botwie2_**: number of threads used by Bowtie2. Default: 6.

### Outputs
- **Kallisto results** (kallisto_ folder) in `Processed_Data/1_Preprocessing/mRNA__kallisto_output`


# Quality Controls

## MRNA_QC__running_fastqc

### Description
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is runned to perform standard quality controls checks on sequenced data.

### Parameters
- **_params.nb_threads_fastqc_**: number of threads used by FastQC. Default: 2.

### Outputs
- **Reads quality control reports** (.zip and .html files) in `Processed_Data/1_Preprocessing/mRNA__fastqc`


## MRNA_QC__running_MultiQC

### Description
A MultiQC html report is made that aggregates all basic FastQC quality controls files.

### Outputs
- **MultiQC report** (.html files) in `Figures_Individual/1_Preprocessing`.
 

