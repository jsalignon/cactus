

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Options](/docs/3_Inputs/Options.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



![](/docs/images/3_Inputs.png "Inputs")

This section covers the inputs files needed by Cactus. There are 3 kind of inputs files: 
 - [Data](/docs/3_Inputs/Data.md) (*.fastq.gz* files): raw sequencing output files
 - [Design](/docs/3_Inputs/Design.md) (*.tsv* files): to indicate the design of the experiment; that is how fastq files relate to samples and conditions, and comparisons to perform and groups of comparisons to plot together in the heatmaps
 - [Options](/docs/3_Inputs/Options.md) (*.config* files): to indicate the parameters to use for the current analysis run. This is the only needed input file for a cactus call.
 
Here is an example of input files at the run directory folder:
<!-- tree -I "results|work"  (worm folder; then editing the output manually)-->
```
.
├── options
│   └── run.yml
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

>**_Note_:** There is an additional .cactus.config file that is located in the root folder and that indicates the global configuration of cactus for all runs of the user.  

>**_Note_:** Directory structure can be changed arbitrarily as files path are specified in the *.yml* input file.

>**_Note_:** The Cactus run will create two additional directories: the results directory and the work directory (a temporary directory created by Nextflow).  
