
* Introduction: [Cactus](./README.md#Cactus), [Quick Start](./README.md#Quick-Start), [Citation](./README.md#Citation)
* Overview: [Short summary](./docs/overview.md#Short-Summary), [Graph](./docs/overview.md#Graph)
* Usage: [Install](./docs/usage.md#Install), [Run](./docs/usage.md#Run), [Reproducibility](./docs/usage.md#Reproducibility), [Troubleshooting](./docs/usage.md#Troubleshooting)
* Input: [Fastq](./docs/input.md#Fastq), [Design](./docs/input.md#Design), [Config](./docs/input.md#Config), [Troubleshooting](./docs/input.md#Troubleshooting)
* Output: [Types](./docs/output.md#Types), [Structure](./docs/output.md#Structure)



# Fastq

 - describe fastq files names


# Design

Inputs of the program are raw ATAC-Seq or mRNA-Seq fastq.gz files. If there are multiple files for the same conditions, the read do not necessarily need to be merged as this can be done within cactus. 
The configuration files are expected by Cactus are: atac_seq.config, mRNA_Seq.config, comparisons.config, regions_to_remove.config, grouped_conditions_for_plots.config, base.config. All are tab separated files. Below are instructions and examples on how to write these files [note: add an example per type of file]:
-	**atac_seq.config**: one line per ATAC-Seq sample. The first entry is the sample ID, and the second is the path were the files is stored. There can be multiple paths, for multiple files, in which cases these files will be merged by Cactus.
-	**mRNA_Seq.config**: one line per mRNA-Seq sample. Same formatting as for atac_seq.config
-	**regions_to_remove.config**: This file allows to remove all reads that map to the specified regions. This is useful in particular for experiments involving RNA interference, as this is known to induce a very strong sequencing signal for the locus that is repressed. The file is configurated this way:
-	**comparisons.config**: this file allows to determine which pairs of conditions will be compared to one another for Differential Abundance Analysis. It is formatted with one comparison per line with one entry per condition. 
-	**grouped_conditions_for_plots.config**: This file allows to define groups of comparisons to plot together on the heatmaps. The format is one line per group with the group ID as first entry and the comparisons as the remaining entries. The comparisons are named this way: comparison1_vs_comparison2.
-	**base.config** (optional): this file allows to overwrite Cactus default settings with custom settings. All parameters from the run.config file can be set changed here to determine how a given experiment is analyzed. See the parameters section below for more details.


# Config

