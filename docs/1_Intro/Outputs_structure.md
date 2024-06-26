
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Tutorial](/docs/1_Intro/tutorial.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Analysis](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment Analysis](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)


<!-- 
homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus

res_dir=$cactus_dir/testing2/singularity/worm/results/

tree -I Run_Info -d -L 3 $cactus_dir/test_datasets/human/results/full_test
tree -I Run_Info -d -L 3 $cactus_dir/testing2/singularity/worm/results/almost_full_test

-->

```
.
├── Figures_Individual
│   ├── 1_Preprocessing
│   │   ├── ATAC__peaks__average_profile
│   │   ├── ATAC__peaks__coverage
│   │   ├── ATAC__peaks__grouped_plots
│   │   ├── ATAC__peaks__saturation_curve
│   │   ├── ATAC__reads__correlations
│   │   ├── ATAC__reads__coverage
│   │   ├── ATAC__reads__insert_size
│   │   └── ATAC__reads__PCA
│   ├── 2_Differential_Analysis
│   │   ├── ATAC__FDR_by_PA
│   │   ├── ATAC__other_plots
│   │   ├── ATAC__PCA_1_2
│   │   ├── ATAC__PCA_3_4
│   │   ├── ATAC__volcano
│   │   ├── mRNA__other_plots
│   │   ├── mRNA__PCA_1_2
│   │   ├── mRNA__PCA_3_4
│   │   ├── mRNA__volcano
│   │   ├── Venn_diagrams__four_ways
│   │   └── Venn_diagrams__two_ways
│   └── 3_Enrichment_Analysis
│       ├── Barplots__CHIP
│       ├── Barplots__chrom_states
│       ├── Barplots__func_anno_BP
│       ├── Barplots__func_anno_KEGG
│       ├── Barplots__genes_self
│       ├── Barplots__motifs
│       ├── Barplots__peaks_self
│       ├── Heatmaps__CHIP
│       ├── Heatmaps__chrom_states
│       ├── Heatmaps__func_anno_BP
│       ├── Heatmaps__func_anno_KEGG
│       ├── Heatmaps__genes_self
│       ├── Heatmaps__motifs
│       └── Heatmaps__peaks_self
├── Figures_Merged
│   ├── 1_Preprocessing
│   ├── 2_Differential_Analysis
│   └── 3_Enrichment_Analysis
├── Processed_Data
│   ├── 1_Preprocessing
│   │   ├── ATAC__peaks__annotated_rds
│   │   ├── ATAC__peaks__split__no_BL_input_RNAi
│   │   ├── ATAC__reads__bam
│   │   ├── ATAC__reads__bam_no_lowQ
│   │   ├── ATAC__reads__bam_no_lowQ_dupli
│   │   ├── ATAC__reads__bam_no_lowQ_dupli_mito
│   │   ├── ATAC__reads__fastq_trimmed
│   │   ├── mRNA__fastqc
│   │   └── mRNA__kallisto_output
│   ├── 2_Differential_Analysis
│   │   ├── ATAC__all_peaks__bed
│   │   ├── ATAC__all_peaks__ChIPseeker
│   │   ├── ATAC__all_peaks__dataframe
│   │   ├── ATAC__all_peaks__DiffBind
│   │   ├── ATAC__all_peaks__gRange
│   │   ├── ATAC__non_annotated_peaks
│   │   ├── DA_split__bed_regions
│   │   ├── DA_split__genes_rds
│   │   ├── mRNA__all_genes__bed_promoters
│   │   └── mRNA__all_genes__rsleuth
│   └── 3_Enrichment_Analysis
│       ├── Barplots__CHIP
│       ├── Barplots__chrom_states
│       ├── Barplots__func_anno_BP
│       ├── Barplots__func_anno_KEGG
│       ├── Barplots__genes_self
│       ├── Barplots__motifs
│       ├── Barplots__peaks_self
│       ├── Enrichment__CHIP
│       ├── Enrichment__chrom_states
│       ├── Enrichment__func_anno_BP
│       ├── Enrichment__func_anno_KEGG
│       ├── Enrichment__genes_self
│       ├── Enrichment__motifs
│       ├── Enrichment__motifs__raw
│       ├── Enrichment__peaks_self
│       ├── Heatmaps__CHIP
│       ├── Heatmaps__chrom_states
│       ├── Heatmaps__func_anno_BP
│       ├── Heatmaps__func_anno_KEGG
│       ├── Heatmaps__genes_self
│       ├── Heatmaps__motifs
│       └── Heatmaps__peaks_self
├── Tables_Individual
│   ├── 1_Preprocessing
│   ├── 2_Differential_Analysis
│   │   ├── ATAC_detailed
│   │   ├── mRNA_detailed
│   │   ├── res_filter
│   │   └── res_simple
│   └── 3_Enrichment_Analysis
│       ├── CHIP
│       ├── chrom_states
│       ├── func_anno_BP
│       ├── func_anno_KEGG
│       ├── genes_self
│       ├── motifs
│       └── peaks_self
└── Tables_Merged
    ├── 1_Preprocessing
    ├── 2_Differential_Analysis
    └── 3_Enrichment_Analysis

```


Cactus produces three main types of outputs: processed data, figures and tables. Both figures and tables are available as merged files (for easier scrolling) or as individual files (for lower memory usage). Tables can be either csv files or formatted Excel files.
