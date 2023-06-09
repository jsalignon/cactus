
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


Cactus was built with the goal to make the installation as pain-free as possible. This was achieved by using only tools within containers or virtual environments for all analysis. Installation of most tools is done automatically the first time the pipeline is run. However, two key dependencies should still be installed: [Nextflow](https://doi.org/10.1038/nbt.3820), and a package manager; one of these tools: [SingularityCE](https://doi.org/10.1371/journal.pone.0177459), [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241), [conda](https://docs.anaconda.com/anaconda/reference/release-notes/) and [Mamba](https://medium.com/@QuantStack/open-software-packaging-for-science-61cecee7fc23). Among these, Singularity is the recommended package manager to use if possible.

The current stable release of Cactus (v0.8.4) was developed and tested using these versions of the dependencies:
 - Nextflow: 22.10.8.5859
 - SingularityCE: 3.11.3-focal
 - Docker: 20.10.21-0ubuntu1~20.04.2
 - conda/Mamba: currently, running Cactus using conda or Mamba is not possible due to a [conda bug](https://github.com/conda/conda/issues/12786). 

In case of issues with running Cactus, please make sure you are using the same version of the dependencies.

You can find here installation links for: [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html#installation), [SingularityCE](https://docs.sylabs.io/guides/latest/admin-guide/installation.html), [Docker](https://docs.docker.com/get-docker/), [conda](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html) and [Mamba](https://mamba.readthedocs.io/en/latest/installation.html).

Please, note that Cactus is currently written in DSL1, which is [not supported anymore](https://nextflow.io/podcast/2023/ep9_end_of_dsl1_chatting_to_bots.html) by the latest default versions of Nextflow. For now, Cactus users should use Nextflow version 22.10.X or earlier, such as v22.10.8 released on April 15 2023. This can be achieved with these commands:
```
wget -qO- https://github.com/nextflow-io/nextflow/releases/download/v22.10.8/nextflow-22.10.8-all | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin
```
