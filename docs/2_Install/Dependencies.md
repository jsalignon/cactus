

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Fastq](/docs/3_Inputs/Fastq.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



Cactus was built with the goal to make the installation as pain-free as possible. This was achieved by using only tools within containers for all analysis. Thereby, installation of most tools is done by simply downloading containers. However, two key dependencies are still necessary. These are [Nextflow](https://doi.org/10.1038/nbt.3820) (the pipeline tool) and [SingularityCE](https://doi.org/10.1371/journal.pone.0177459) (the container tool).  

To install these tools please follow the instructions on these links: [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html), [SingularityCE](https://docs.sylabs.io/guides/latest/admin-guide/installation.html).  

In case of issue with the pipeline and/or for reproducibility purposes, it might be useful to try to run Cactus with the version of these tools that was used when developing the current Cactus release. These are:  
  - SingularityCE: version 3.10.0+91-g13f189977 (released on May 17, 2022)  
  - Nextflow: version 22.05.0-edge build 5704 (released on May 25, 2022 )  
