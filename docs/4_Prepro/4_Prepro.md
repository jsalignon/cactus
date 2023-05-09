
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



This section covers standard preprocessing analysis. This corresponds to peak calling for ATAC-Seq and to transcripts abundance quantification for mRNA-Seq, as well as quality control analysis.

The pipeline was originally developed to analyze ATAC-Seq data and therefore more quality controls are being performed for ATAC-Seq data. An important difference between ATAC-Seq and mRNA-Seq data preprocessing is that reads are aligned for ATAC-Seq data analysis, while an alignment-free method is used for mRNA-Seq. For this reason the preprocessing is much faster for mRNA-Seq data. This also means that less quality controls can be made for mRNA-Seq data. More quality controls options for mRNA-Seq data may be added if needed in future versions of Cactus.

For mRNA-Seq data, transcripts abundances are quantifyied with [kallisto](https://doi.org/10.1038/nbt.3519) and quality controls are made with [FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) and [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/).

ATAC-Seq preprocessing steps mostly follow the [guidelines (first version)](https://informatics.fas.harvard.edu/atac-seq-guidelines-old-version.html) from the Harvard Faculty of Arts and Sciences. With the key steps being: reads merging, trimming, aligning, filtering (low quality, duplicates, mitonchondrial and small contigs), shifting (transposase-shift), and peaks calling (with [MACS2](https://doi.org/10.1186/gb-2008-9-9-r137)), splitting, and filtering (blacklisted regions, input control, specific regions). Quality controls are made via published tools ([FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/), [MultiQC](https://pubmed.ncbi.nlm.nih.gov/27312411/), [DeepTools](https://doi.org/10.1093/nar/gkw257) for reads profiles and correlation, [ChIPseeker](http://dx.doi.org/10.1093/bioinformatics/btv145) for distribution of annotated peaks) and homemade scripts (saturation curve, reads overlap with genomic regions, ...).
