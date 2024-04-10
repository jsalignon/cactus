

# Menu
- [Configuration files](#Configuration-files)
- [Mandatory parameters](#Mandatory-parameters)

# Introduction

This tutorial will highlight some use case of Cactus using the worm test dataset.
The full script for this tutoral can be found here[TO DO].
We will use Singularity for this tutorial (see the [dependencies section](/docs/2_Install/Dependencies.md) if needed).

# Setting up global parameters, downloading references, and test datasets

We should first create a Cactus configuration file. For instance, this can be done using this command (we disable [Nextflow Tower](https://cloud.tower.nf/) for this example):
```bash
cat >> ${HOME}/.cactus.config << EOL
params.references_dir        = "${HOME}/workspace/cactus/references"
params.singularity_cache_dir = "${HOME}/workspace/singularity_containers"
params.tower_token           = "*"
params.enable_tower          = false
EOL
```

We can export the nextflow version to make sure we use the correct:
```bash
export NXF_VER=22.10.8
```

For reproducibility purposes, it is useful to specify a version of Cactus whenever running an analysis. We will use the latest version (v0.9.0) for this tutorial.
```bash
cactus_version=0.9.0
```

We can then download the references and test datasets using this command:
```bash
nextflow run jsalignon/cactus/scripts/download/download.nf -r $cactus_version --test_datasets --references -profile singularity --species worm
```

This command will download the references to the `params.references_dir` folder set up in the global Cactus configuration file, as well as the worm test dataset in the current folder.

The worm test dataset comes from [TODO].

We can go to the worm test dataset folder and check the pre-made experimental design with these commands:
```bash
cd worm
printf "\n\nATAC-Seq fastq files\n"
cat design/atac_fastq.tsv
printf "\n\nmRNA-Seq fastq files\n"
cat design/mrna_fastq.tsv
printf "\n\nComparisons to be made\n"
cat design/comparisons.tsv
printf "\n\nGroups of comparisons to be plot together in heatmaps\n"
cat design/groups.tsv
```

Which gives:
```
ATAC-Seq fastq files
input   data/atac/sample_200K_reads_atac_SRX2333004_SRR5000684_R1.fastq.gz
ctl_1   data/atac/sample_200K_reads_atac_SRX3029124_SRR5860424_R1.fastq.gz
ctl_2   data/atac/sample_200K_reads_atac_SRX3029125_SRR5860425_R1.fastq.gz
ctl_3   data/atac/sample_200K_reads_atac_SRX3029126_SRR5860426_R1.fastq.gz
spt16_1 data/atac/sample_200K_reads_atac_SRX3029130_SRR5860430_R1.fastq.gz
spt16_2 data/atac/sample_200K_reads_atac_SRX3029131_SRR5860431_R1.fastq.gz
spt16_3 data/atac/sample_200K_reads_atac_SRX3029132_SRR5860432_R1.fastq.gz
hmg4_1  data/atac/sample_200K_reads_atac_SRX3029133_SRR5860433_R1.fastq.gz
hmg4_2  data/atac/sample_200K_reads_atac_SRX3029134_SRR5860434_R1.fastq.gz
hmg4_3  data/atac/sample_200K_reads_atac_SRX3029135_SRR5860435_R1.fastq.gz

mRNA-Seq fastq files
ctl_1   data/mrna/sample_50K_reads_mrna_SRX3029112_SRR5860412.fastq.gz
ctl_2   data/mrna/sample_50K_reads_mrna_SRX3029113_SRR5860413.fastq.gz
ctl_3   data/mrna/sample_50K_reads_mrna_SRX3029114_SRR5860414.fastq.gz
hmg4_1  data/mrna/sample_50K_reads_mrna_SRX3029115_SRR5860415.fastq.gz
hmg4_2  data/mrna/sample_50K_reads_mrna_SRX3029116_SRR5860416.fastq.gz
hmg4_3  data/mrna/sample_50K_reads_mrna_SRX3029117_SRR5860417.fastq.gz
spt16_1 data/mrna/sample_50K_reads_mrna_SRX3029118_SRR5860418.fastq.gz
spt16_2 data/mrna/sample_50K_reads_mrna_SRX3029119_SRR5860419.fastq.gz
spt16_3 data/mrna/sample_50K_reads_mrna_SRX3029120_SRR5860420.fastq.gz

Comparisons to be made
hmg4 ctl
spt16 ctl
hmg4 spt16

Groups of comparisons to be plot together in heatmaps
all     hmg4_vs_ctl     spt16_vs_ctl    hmg4_vs_spt16
ctl     hmg4_vs_ctl     spt16_vs_ctl
spt16   spt16_vs_ctl    hmg4_vs_spt16
```

We will not remove any genes and regions at first, and use the empty files `design/regions_to_remove_empty.tsv` and `design/genes_to_remove_empty.tsv`.

We can now run Cactus. This can be done using the parameter file, such as the `parameters/full_test.yml` file, which contains:
```
res_dir                   : 'results/full_test'
species                    : 'worm'
chromatin_state           : 'iHMM.M1K16.worm_L3'
split__threshold_type     : 'rank'
split__threshold_values   : [ 200, 1000 ]
```
However, for this tutorial, we will set up the parameters directly in the command line for better clarity. Nextflow parameters are set up with a single hyphen (i.e., "-profile"), while Cactus parameters are set up with double hyphens (i.e., "--res_dir").

It is a good idea to run Cactus sequentially, to make sure the different steps of analysis are correct, and to avoid the need to re-run the whole pipeline again. Since Cactus is built with the Nextflow language, caching is efficiently implemented. Therefore, changing one paramter or one sample will only affect the specific process or sample that was changed. 

We therefore run Cactus by using the `--params.disable_all_enrichments` parameter to stop the analysis after the second step (differential analysis) with the following command:
```bash
nextflow run jsalignon/cactus -profile singularity -r $cactus_version \
	--res_dir results/tutorial              \ # the output folder
	--species worms                         \ # species under study (mandatory)
	--chromatin_state iHMM.M1K16.worm_L3    \ # chromatin state file to use (mandatory)
	--split__threshold_type rank            \ # using rank(-log10(FDR)) cuttoff
	--split__threshold_values [ 200, 1000 ] \ # selecting the top 200 or 1000 most significant results
	--disable_all_enrichments               \ # disabling all enrichment analysis

```

Note, that creation of bigWig files, as well as the saturation curve analysis can both be computationally intensive and therefore there are parameters to disable these steps if needed (`params.do_bigwig` and `params.do_saturation_curve`).

After this step, we can look at the results from the preprocessing and quality control analyses to ensure all samples are of good enough quality. This can be assessed for both ATAC-Seq and mRNA-Seq samples by looking at the MultiQC report files, present in the folder `results/tutorial/Figures_Merged/1_Preprocessing`. Many additional quality control figures are available for ATAC-Seq in the same folder, including reads and peaks coverage, reads insert size, PCA and correlation matrix. The saturation curve analysis helps to determine if one should resequence or not. This can be assessed by looking at the shape of the curve. For instance, if the curve is still steeply increasing, then it might be needed to sequence the samples further.




-params-file parameters/full_test.yml -r main -latest

res_dir                   : 'results/full_test'
species                    : 'worm'
chromatin_state           : 'iHMM.M1K16.worm_L3'
split__threshold_type     : 'rank'
split__threshold_values   : [ 200, 1000 ]



```
cat design/regions_to_remove_empty.tsv
```

We will now run Cactus. We don't specify the peaks to remove at first.


# Setting up configuration files

The previous

# Running preprocessing and checking quality controls

# Running differential analysis, filtering, and selecting PA filters

# Running enrichment analysis between comparisons 

# Running the full analysis
