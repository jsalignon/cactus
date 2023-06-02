
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


Cactus was built with the goal to make the installation as pain-free as possible. This was achieved by using only tools within containers or virtual environments for all analysis. Installation of most tools is done automatically the first time the pipeline is run. However, two key dependencies should still be installed. These are:

1. The workflow language [Nextflow](https://doi.org/10.1038/nbt.3820) ([Install](https://www.nextflow.io/docs/latest/getstarted.html#installation)).
Note that Cactus is currently written in DSL1, which is not supported by the latest default versions of Nextflow. Cactus users should use Nextflow version 22.10.X or earlier, such as v22.10.8 released on April 15 2023. This can be achieved with these commands:

```
wget -qO- https://github.com/nextflow-io/nextflow/releases/download/v22.10.8/nextflow-22.10.8-all | bash
chmod +x nextflow
sudo mv nextflow /usr/local/bin
```

2. A tools manager that can be any of:
 - [SingularityCE](https://doi.org/10.1371/journal.pone.0177459) ([Install](https://docs.sylabs.io/guides/latest/admin-guide/installation.html))
 - [Docker](https://dl.acm.org/doi/10.5555/2600239.2600241) ([Install](https://docs.docker.com/get-docker/))
 - [conda](https://docs.anaconda.com/anaconda/reference/release-notes/) ([Install](https://docs.conda.io/projects/conda/en/latest/user-guide/install/index.html))
 - [Mamba](https://medium.com/@QuantStack/open-software-packaging-for-science-61cecee7fc23) ([Install](https://mamba.readthedocs.io/en/latest/installation.html))

When using conda or mamba, the set-up should be as described in [Bioconda's usage section](https://bioconda.github.io/#usage):
```
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

During Cactus development, Mamba was installed with this command: `conda install mamba -n base -c conda-forge`

In case of issue with the pipeline and/or for reproducibility purposes, it might be useful to try to run Cactus with the version of these tools that were used when developing the current Cactus release. These are:  
  - Nextflow: 22.09.2-edge.5765
  
  - SingularityCE: 3.10.0+91-g13f189977
  - Conda: 22.9.0
  - Mamba 0.27.0 (conda 22.9.0)

Docker version 18.06.1-ce, build e68fc7a
