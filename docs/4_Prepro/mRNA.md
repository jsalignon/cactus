
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

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

### Parameters
- **_params.kallisto__nb_threads_**: number of threads used by kallisto. Default: 6.
- **_params.kallisto__bootstrap_**: Number of bootstrap samples. Default: '100'.
- **_params.kallisto__fragment_len_**: Estimated average fragment length. For single end only. Default: '180'.
- **_params.kallisto__fragment_sd_**: Estimated standard deviation of fragment length. For single end only. Default: '20'.

### Outputs
- **Kallisto results**: `Processed_Data/1_Preprocessing/mRNA__kallisto_output/kallisto_${sample}/{abundance.{h5,tsv},run_info.json}`.


# Quality Controls

## MRNA_QC__running_fastqc

### Description
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) is runned to perform standard quality controls checks on sequenced data.

### Parameters
- **_params.fastqc__nb_threads_**: number of threads used by FastQC. Default: 2.

### Outputs
- **Reads quality control reports**: `Processed_Data/1_Preprocessing/mRNA__fastqc/*_R{1,2}_fastqc.html`.




## MRNA_QC__running_MultiQC

### Description
A MultiQC html report is made that aggregates all basic FastQC quality controls files.

### Outputs
- **MultiQC report**: `mRNA__multiQC.html`. [Example](https://htmlpreview.github.io/?https://github.com/jsalignon/cactus/blob/main/docs/examples/html/mRNA__multiQC.html).
 
### Output folders
- `Figures_Individual/1_Preprocessing`
- `Figures_Merged/1_Preprocessing`
