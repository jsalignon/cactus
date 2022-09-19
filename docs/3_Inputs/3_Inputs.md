

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [fastq](/docs/3_Inputs/fastq.md), [tsv](/docs/3_Inputs/tsv.md), [config](/docs/3_Inputs/config.md), [yml](/docs/3_Inputs/yml.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



![](/docs/images/3_Inputs.png "Inputs")

This section covers the inputs files needed by Cactus. There are 4 kind of inputs files: 
 - [data](/docs/3_Inputs/Inputs_data.md): raw *.fastq.gz* sequencing output files
 - [design](/docs/3_Inputs/Design.md): a set of *.tsv* files to indicates how fastq files relate to samples, conditions, comparisons to perform and groups of comparisons to plot together in heatmaps
 - [configuration](/docs/3_Inputs/Configuration.md): two *.config* files that indicates the custom parameters of cactus to use for the current analysis run
 - [configuration](/docs/3_Inputs/Configuration.md): the *.yml* input file to a cactus call
 
Here is an overview of the input files at the run directory folder:
<!-- tree -I "results|work"  (worm folder; then editing the output manually)-->
```
.
├── conf
│   └── run.config
├── data
│   ├── atac
│   │   ├── sample_1000K_reads_atac_SRX2333004_SRR5000684_R1.fastq.gz
│   │   ├── sample_1000K_reads_atac_SRX2333004_SRR5000684_R2.fastq.gz
│   │   ├── sample_1000K_reads_atac_SRX3029124_SRR5860424_R1.fastq.gz
│   │   ├── ...
│   └── mrna
│       ├── sample_50K_reads_mrna_SRX3029112_SRR5860412.fastq.gz
│       ├── sample_50K_reads_mrna_SRX3029113_SRR5860413.fastq.gz
│       ├── ...
└── design
    ├── atac_fastq.tsv
    ├── comparisons.tsv
    ├── groups.tsv
    ├── mrna_fastq.tsv
    └── regions_to_remove.tsv
```

>**_Note_:** There is an additional and optional global.config file that is located in the cactus folder and indicate the global configuration of cactus for all runs of the user.

>**_Note_:** Input data files do not necessarily be located at the same place as the analysis directory as shown in this directory tree.

>**_Note_:** The Cactus run will create two additional directories: the results directory and the work directory (a temporary directory created by Nextflow)
