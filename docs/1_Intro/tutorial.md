
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



# Menu
- [Introduction](#Introduction)
- [Setting up global parameters, downloading references, and test datasets](#Setting-up-global-parameters,-downloading-references,-and-test-datasets)
- [Running Cactus without enrichment analysis](#Running-Cactus-without-enrichment-analysis)
- [Checking quality controls and removing artifical signal](#Checking-quality-controls-and-removing-artifical-signal)
- [Selecting Differential Analysis Subset filters](#Selecting-Differential-Analysis-Subset-filters)
- [Looking at correlation between comparisons for different experiment types](#Looking-at-correlation-between-comparisons-for-different-experiment-types)
- [Running the full analysis](#Running-the-full-analysis)


# Introduction

This tutorial will highlight some use case of Cactus using the worm test dataset.
The full script for this tutoral can be found [here](/docs/1_Intro/tutorial.sh).

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

We can export the nextflow version to make sure we use a compatible version:
```bash
export NXF_VER=22.10.8
```

For reproducibility purposes, it is useful to specify a version of Cactus whenever running an analysis. We will use the latest available version (v0.9.0) for this tutorial:
```bash
cactus_version=0.9.0
```

We can then download the references and test datasets using this command:
```bash
nextflow run jsalignon/cactus/scripts/download/download.nf -r $cactus_version --test_datasets --references -profile singularity --species worm
```

This command will download the references to the `params.references_dir` folder set up in the global Cactus configuration file, as well as the worm test dataset in the current folder.

The worm test dataset comes from the Koluntzic et al manuscript (see the [Test datasets](/docs/2_Install/Test_datasets.md) section), and contains ATAC-Seq and RNA-Seq samples of *Caenorhabditis elegans* worms exposed to RNA interference (RNAi) on two FACT subunits; hmg-4, spt-16, as well as Renilla luciferase (Rluc) which served as a control (abbreviated "ctl").

We can go to the freshly downloaded worm test dataset folder and check the pre-made experimental design with these commands:
```bash
cd worm
printf "\n\nATAC-Seq fastq files\n"
cat design/atac_fastq.tsv
printf "\n\nmRNA-Seq fastq files\n"
cat design/mrna_fastq.tsv
printf "\n\nComparisons to be made\n"
cat design/comparisons.tsv
printf "\n\nGroups of comparisons to be plotted together in heatmaps\n"
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

Finally, we should set up two global parameters for the enrichment analysis part: 

 - The mandatory `--params.chromatin_state` parameter defines the file that will be used to do enrichment of chromatin states. The options can be found in the file `$reference_dir_folder/encode_chromatin_states_metadata.csv` for ChromHMM-derived files and on this [website](https://compbio.med.harvard.edu/modencode/webpage/hihmm/) (or using the command `ls $reference_dir_folder/chromatin_states/`) for HiHMM-derived states. See the [References section](https://github.com/jsalignon/cactus/blob/main/docs/2_Install/References.md#Pipeline-to-get-references) for details on the origin of these files.

 - The optional `--chip_ontology` parameter define the ontology term that will be used to filter the ENCODE ChIP-Seq profiles on which to conduct enrichment analysis. The options can be found in the [References section](https://github.com/jsalignon/cactus/blob/main/docs/2_Install/References.md#CHIP-Seq) (or using the command  `cat $reference_dir_folder/chip_ontology_groups_sizes.txt`).
Information about the details of the ChIP-Seq profiles abbreviations can be found in the file `$reference_dir_folder/encode_chip_metadata.csv`.



# Running Cactus without enrichment analysis

We can now run Cactus. This can be done for instance by adding this argument `-params-file parameters/full_test.yml` to the Cactus call, which will set all of the parameters indicated in the provided YML file. Here is the content of the `parameters/full_test.yml` file:
```
res_dir                   : 'results/full_test'
species                   : 'worm'
chromatin_state           : 'iHMM.M1K16.worm_L3'
split__threshold_type     : 'rank'
split__threshold_values   : [ 200, 1000 ]
```

For this tutorial, we will set up all the parameters directly in the command line. Nextflow parameters are set up with a single hyphen (i.e., "-profile"), while Cactus parameters are set up with double hyphens (i.e., "--res_dir").


It is a good idea to run Cactus sequentially, to make sure the different steps of analysis are correct, and to avoid the need to re-run the whole pipeline again. Since Cactus is built with the Nextflow language, caching is efficiently implemented. Therefore, changing one parameter or one sample will only affect the parts of the analysis that are affected by the change.

We therefore run Cactus and stop the analysis after the second step (differential analysis) with the following command:
```bash
nextflow run jsalignon/cactus -profile singularity -r $cactus_version \
	--res_dir results/tutorial              \
	--species worm                          \
	--chromatin_state iHMM.M1K16.worm_L3    \
	--design__regions_to_remove design/regions_to_remove_empty.tsv \
	--disable_all_enrichments               
```
With: 
 - `--res_dir`: the output folder
 - `--species`: species under study (mandatory)
 - `--chromatin_state`: chromatin state file to use for the chromatin enrichment analysis (mandatory)
 - `--design__regions_to_remove`: using an an empty file that do not contain any regions to remove from ATAC-Seq analyses (more on that below).
 - `--disable_all_enrichments`: disabling all enrichment analyses

Note also, that creation of bigWig files, as well as the saturation curve analysis can both be computationally intensive and therefore there are parameters to disable these steps if needed (`params.do_bigwig` and `params.do_saturation_curve`).


# Checking quality controls and removing artifical signals

After this step, we can look at the results from the preprocessing and quality control analyses to ensure all samples are of good enough quality. This can be assessed for both ATAC-Seq and mRNA-Seq samples by looking at the MultiQC report files, present in the folder `results/tutorial/Figures_Merged/1_Preprocessing`. Many additional quality control figures are available for ATAC-Seq in the same folder, including reads and peaks coverage, reads insert size, PCA and correlation matrix. 

The saturation curve analysis helps to determine if one should resequence or not. This can be assessed by looking at the shape of the curve. For instance, if the curve is still steeply increasing, then it would likely be useful to sequence the samples further. Here is the saturation curve for an hmg-4 RNAi sample:

<img src="/docs/examples/png/hmg4_3__saturation_curve.png" width="300" />

Since this a test datasets with few randomly sampled reads, 
We can see that the saturation curve is still steeply increasing (even though we notice a strange break in the curve). Since a clear plateau has not been reached, this figure indicates that sequencing deeper could be helpful to identify more ATAC-Seq peaks. This is expected as this test datasets contains few randomly sampled reads.

Looking at the ATAC-Seq volcano plot for the hgm-4 vs stp-16 comparison, we can notice highly significant peaks nearby the spt-16 and hmg-4 genes:

<img src="/docs/examples/png/hmg4_vs_spt16__ATAC_volcano__no_rtr.png" width="300" /> 

This is due to RNA interference which induce a very strong sequencing signal at the repressed locus.

Therefore, we repeat our analysis by excluding the regions of the targeted gene for ATAC-Seq samples, using the file `design/regions_to_remove.tsv` which contains:
```
hmg4    Hmg4->chrIII:7,379,143-7,381,596
spt16   Spt16->chrI:10,789,130-10,793,152
```

The new commands is:
```bash
nextflow run jsalignon/cactus -profile singularity -r $cactus_version \
	--res_dir results/tutorial              \
	--species worms                         \
	--chromatin_state iHMM.M1K16.worm_L3    \
	--design__regions_to_remove design/regions_to_remove.tsv \
	--disable_all_enrichments               
```
Looking at the volcano plot, we can see that articial signals at the hmg-4 and spt-16 locuses have been removed:

<img src="/docs/examples/png/hmg4_vs_spt16__ATAC_volcano.png" width="300" />      


# Selecting Differential Analysis Subset filters

At this stage, we can also look at which DAS (Differential Analysis Subset) filters could make sense to use for splitting the Differential Analysis (DA) results into subsets. For instance, we can look at the volcano plots for ATAC-Seq (left image) and mRNA-Seq (right image:

<img src="/docs/examples/png/hmg4_vs_ctl__ATAC_volcano.png" width="300" /> <img src="/docs/examples/png/hmg4_vs_ctl__mRNA_volcano.png" width="300" />

We can see that the p-values are much smaller for mRNA-Seq than for ATAC-Seq. Indeed, only two peaks are significant opened, while many genes are differentially expressed. 

From these plots, one could run the analysis by using different -log10(FDR) cutoffs using these parameters: `--split__threshold_type FDR` and `--split__threshold_values [ 1.3, 3, 10 ]` (please observe that -log10(1.3) ~= 0.05). This would allow to dissect the enrichment of differentially expressed genes at various significance cutoffs. This can be helpful as different strengh of induction of Transcription Factors can lead to different significance of association of their target genes. 

Alternatively, one could try to check if there is consistency between the ATAC-Seq and mRNA-Seq results by selecting the top N most significant DAS, using these parameters: `--split__threshold_type rank` and `--split__threshold_values [ N ]`. 

One can also check if certain Peak Annotation (PA) filters would be interesting to investigate by looking at the FDR by PA plot:

<img src="/docs/examples/png/hmg4_vs_ctl__ATAC_FDR_by_PA.png" width="300" />

In this case, no PA filter seem to be especially enriched in significant peaks, therefore we can keep the default parameter: `--split__peak_assignment ['all']`.

After having identified the PA filters to use, we can launch a new analysis run that conducts enrichment analysis this time. This step can be done very quickly when excluding functional enrichment and motifs enrichment analysis (i.e., by default, all enrichment analysis steps are conducted), as these are more computationally intensive. Indeed, enrichment of genomic regions, like ChIP-Seq, chromatin states, or DASs can be done very rapidly as it is based on overlaping bed files using bedtools. We can now use this command to do a quick enrichment analysis:
```bash
nextflow run jsalignon/cactus -profile singularity -r $cactus_version \
	--res_dir results/tutorial              \
	--species worms                         \
	--chromatin_state iHMM.M1K16.worm_L3    \
	--design__regions_to_remove design/regions_to_remove.tsv \
	--split__threshold_type     rank        \
	--split__threshold_values [ 200, 1000 ] \
	--do_func_anno_enrichment   false       \
	--do_motif_enrichment       false
```
With: 
 - `--split__threshold_type`: using rank(-log10(FDR)) cutoff instead of the default -log10(FDR) cutoff
 - `--split__threshold_values`: selecting the top 200 and 1000 most significant results
 - `--do_func_anno_enrichment`: disabling functional annotation enrichment analysis.
 - `--do_motif_enrichment`:     disabling motifs enrichment analysis.


# Looking at correlation between comparisons for different experiment types

At this point, a key result to look at is the `Heatmaps__genes_self.pdf` and `Heatmaps__peaks_self.pdf` figures in the `results/tutorial/Figures_Merged/3_Enrichment_Analysis` folder. These figures will give a concise overview of how the various comparisons correlate with each other. Here are the gene self-overlap figures for the top 1000 most significant Differential Analysis (DA) entries per comparison, with the DEGs (Differentially Expressed Genes) on the left, the genes associated with DARs (Differentially Accessible regions) in the middle, and the HA-HE and LA-LE genes (DEGs with neighbouring DARs with expected chromatin accessibility changes) on the right:

<img src="/docs/examples/png/ATAC__all__1000__ctl__genes_self__heatmap.png" width="233" /> <img src="/docs/examples/png/mRNA__Null__1000__ctl__genes_self__heatmap.png" width="233" /> <img src="/docs/examples/png/both__all__1000__ctl__genes_self__heatmap.png" width="233" />

These figures reveal a clear clustering between changes induced by the hmg-4 and the spt-16 knockdown in both experiment type, which is consistent with them being member of the same complex (the FACT complex). 

The similar number of DA entries is expected as we selected the top 1000 most significant DA entries per comparisons. However, it is striking to observe the same number of shared repressed regions/genes between knockdown of the hmg-4 and spt-16 in the ATAC-Seq and mRNA-Seq analysis (~290 each), while the ATAC-Seq pvalues were much higher. 

This result support the idea that p-values can not always directly be compared between assays and experiment to infer true functional/biological relevance of the results. And therefore, it can make sense to focus on the top N most significant entries between the two omics experimental results to study consistent chromatin and gene expression changes. 
Interestingly, while only ~30 repressed genes had a nearby closing in chromatin for both factors, we can see that a fifth of those are overlaping (adjusted p-value of 0.000015). 

Altogether, these results indicate that using p-vaule rank cutoffs can be a useful way to combine the two omics data type in the case of widely different significance levels. However, sequencing deeper the ATAC-Seq data would obviously be the recommended approach to have more robust ATAC-Seq results.


# Running the full analysis

Finally, we can complete the full analysis with this command:
```bash
nextflow run jsalignon/cactus -profile singularity -r $cactus_version \
	--res_dir results/tutorial              \
	--species worms                         \
	--chromatin_state iHMM.M1K16.worm_L3    \
	--design__regions_to_remove design/regions_to_remove.tsv \
	--split__threshold_type     rank        \
	--split__threshold_values [ 200, 1000 ]
```

Since previous analyses are cached, this will only run the analysis steps that were not made before, which are functional annotations (defaults: KEGG and GO-BP) and motifs enrichment analyses.

Here are the KEGG enrichment results for ATAC-Seq (left), mRNA-Seq (middle) and both (right) from the `results/tutorial/Figures_Merged/3_Enrichment_Analysis` folder:

<img src="/docs/examples/png/ATAC__all__1000__ctl__func_anno_KEGG__heatmap.png" width="233" /> <img src="/docs/examples/png/mRNA__Null__1000__ctl__func_anno_KEGG__heatmap.png" width="233" /> <img src="/docs/examples/png/both__all__1000__ctl__func_anno_KEGG__heatmap.png" width="233" />

We can see that genes associated with DARs are enriched in few pathways. 

At the transcriptional level, we can see a large up-regulation of various metabolism pathways upon loss of FACT subunits.

Finally, High Accessibility-High Expression (HA-HE) genes upon knockdown of spt-16 are also enriched in metabolism pathways, while Low Accessibility-Low Expression (LA-LE) genes upon knockdown of both FACT subunit are enriched in protein export and phagosome.
Of course, the biological interpretation of these results are very limited as the test datasets analysed contains few reads.

To have more details on the enrichment results, one can look at the detailed Excel (or csv) tables in the `results/tutorial/Tables_Merged/3_Enrichment_Analysis` folder. Filtering the KEGG table for `Experiment Type (ET) = "both"` (which keeps only HA-HE and LA-LE gene sets), we obtain:

<img src="/docs/examples/xlsx_png/both__all__1000__ctl__KEGG__table.png" width="700" />

In this table, we can see that 3 (`ov_da` column) out of 11 (`tot_da` column) LA-LE genes upon hmg-4 knockdown are enriched in the "Phagosome" KEGG pathway (and these genes are indicated in the `genes_id` column). Therefore ~30 percent of HA-HE genes (`pt_da` column) are members of this pathway, in contrast to only 2% of the background genes (`pt_nda` column). This result in a p-value of 0.0015 (`pval` column), an adjusted p-value of 0.026 (`padj` column), and a Log2 Odds Ratio of 4.11 (`L2OR` column). Please note that the table indicates only 11 LA-LE genes for the hmg-4 vs control comparison, while there are in total 29 such genes, as shown in the heatmap of the correlation between comparisons above. This is because only genes that are present in the KEGG database are being considered for KEGG enrichment analyses.

***Additional observations:***
 
 - Please, note that the pipeline can currently break if the heatmap filtering parameters reduce certain DAS to 1 entry or less. If such crash occur, one can modify the `--heatmaps_filter__func_anno`, `--heatmaps_filter__CHIP`, and the `--heatmaps_filter__motifs` filters to try to keep more entries and avoid the crash.

 - The `--barplots_ggplot__*`, `--heatmaps_params__*`, and `--heatmaps_ggplot__*` parameters can be used to adjust the selection of terms to display as well as the general appearance of the figures (see the [parameters section](/docs/3_Inputs/Parameters.md) for details). 
