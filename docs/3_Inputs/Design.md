

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [Data](/docs/2_Install/Install_data.md), [Test_datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Inputs_data.md), [Design](/docs/3_Inputs/Design.md), [Configuration](/docs/3_Inputs/Configuration.md)
* [Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [Input Files](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [Differential Abundance](/docs/5_DA/5_DA.md): [DBA](/docs/5_DA/DBA.md), [DGEA](/docs/5_DA/DGEA.md), [Split](/docs/5_DA/Split.md), [Outputs](/docs/5_DA/Outputs.md)
* [Enrichment](/docs/6_Enrich/6_Enrich.md): [Overlap](/docs/6_Enrich/Overlap.md), [Outputs](/docs/6_Enrich/Outputs.md)

[](END_OF_MENU)



### atac_fastq.tsv: fastq files for ATAC-Seq
**Format:** *sample_id fastq_file_path*
**Description:** Each line referes to one fastq file. If two files have the same sample_id they will be considered to be different sequencing runs of the same sample and will be merged by Cactus. 
**Fields details:**
 - *sample_id*: The sample_id is a combination of a condition_id and a replicate_number, united with an underscore. Note that condition_ids can only contain alphanumerical characters (A-Z, a-z and 0-9), and cannot contain special characters (such as underscore). 
 - *fastq_file_path*: can be either an absolute path or a relative path (recommended) from the directory where Cactus is run

**Example:**
```
ctl_1 data/atac/sample_200K_reads_atac_SRX3029124_SRR5860424_R1.fastq.gz
ctl_2 data/atac/sample_200K_reads_atac_SRX3029125_SRR5860425_R1.fastq.gz
ctl_2 data/atac/sample_200K_reads_atac_SRX3029125_SRR5860426_R1.fastq.gz
hmg4_1 data/atac/sample_200K_reads_atac_SRX3029133_SRR5860433_R1.fastq.gz
hmg4_2 data/atac/sample_200K_reads_atac_SRX3029134_SRR5860434_R1.fastq.gz
spt16_1 data/atac/sample_200K_reads_atac_SRX3029130_SRR5860430_R1.fastq.gz
spt16_2 data/atac/sample_200K_reads_atac_SRX3029131_SRR5860431_R1.fastq.gz
```


### mrna_fastq.tsv: fastq files for mRNA-Seq
Format: Same formatting as for atac_fastq.tsv

### regions_to_remove.tsv: regions to filter out during ATAC-Seq peaks preprocessing
Format: *condition_id, Locus_name->genomic_coordinates (chromosome:start-end)*

Peaks that overlap with regions defined in this file will be excluded from the analysis. This is particularly useful in experiments involving RNA interference, as this is known to induce a very strong sequencing signal at the repressed locus. The first field indicates in which sample the region should be removed and the second is the coordinates of the region to remove. Multiple regions can be removed for the same sample by adding multiple lines.

Example:
```
gaf gaf->3L:14,747,929-14,761,049
b170 bap170->2R:6,636,512-6,642,358
n301 nurf301->3L:233,926-246,912
n301b170 bap170->2R:6,636,512-6,642,358
n301b170 nurf301->3L:233,926-246,912
```

**__Note__:** this file can be empty if no region needs to be removed

### comparisons.tsv: pairs of conditions to compare during Differential Abundance Analysis
This file is formatted with one comparison per line with one entry per condition. 
second sample is the "control"
Example:
```
hmg4 ctl
spt16 ctl
hmg4 spt16
```

### groups.tsv: groups of comparisons to plot together on the heatmaps plots. 
The format is one line per group with the group ID as first entry and the comparisons as the remaining entries. The comparisons are named this way: condition1_vs_condition2. 

Example:
```
all hmg4_vs_ctl spt16_vs_ctl hmg4_vs_spt16
ctl hmg4_vs_ctl spt16_vs_ctl
spt16 spt16_vs_ctl hmg4_vs_spt16
```


**__Final Note__:** Fields in tsv design files should be separated by space or tabs.
