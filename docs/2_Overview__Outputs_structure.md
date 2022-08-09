
* Introduction: [Cactus](/README.md), [Quick Start](/docs/1_Intro__Quick_start.md), 
* Overview: [Graph](/docs/2_Overview__Graph.md), [Outputs structure](/docs/2_Overview__Outputs_structure.md)
* Install: [Dependencies](/docs/3_Install__Dependencies.md), [Containers](/docs/3_Install__Containers.md), [Data](/docs/3_Install__Data.md), [Test_datasets](/docs/3_Install__Test_datasets.md)
* Run: [Input Data](/docs/4_Run__Input_data.md), [Input Files](/docs/4_Run__Input_files.md), [Parameters](/docs/4_Run__Parameters.md)
* Preprocessing: ATAC: [Method](/docs/5_AP__Method.md), [Figures](/docs/5_AP__Figures.md), [MultiQC](/docs/5_AP__MultiQC.md), mRNA: [Method](/docs/6_MP__Method.md), [MultiQC](/docs/6_MP__MultiQC.md)
* Differential Abundance: [ATAC](/docs/7_DA__DiffBind.md), [mRNA](/docs/7_DA__Sleuth.md), [Figures](/docs/7_DA__Figures.md), [Tables](/docs/7_DA__Tables.md)
* Splitting peak sets: [Split](/docs/8_SP__Splitting_peak_sets.md), [Venn Diagrams](/docs/8_SP__Venn_diagrams.md)
* Enrichment: [ChIP](/docs/9_Enrich__CHIP.md), [Chromatin state](/docs/9_Enrich__Chromatin_states.md), [Motifs](/docs/9_Enrich__Motifs.md), [Func. Anno.](/docs/9_Enrich__Functional_annotations.md), [Figures](/docs/9_Enrich__Figures.md), [Tables](/docs/9_Enrich__Tables.md)


[](END_OF_MENU)







# Types

[Note: There should be one example Figure for each output]
There are three main types of outputs: processed data, figures and tables. Both figures and tables are available either merged (for easier scrolling) and separated (for lower memory usage). Tables can be either csv files or formatted Excel files.
Here is the structure of the output file with a short description of each output file and examples when relevant: 
TO DO. [add the tree of output as clickable and with a submenu, and description and image for each item]


# Structure

<!-- Tree structure obtained with this command in the results folder: tree -I Run_Info -d -L 3 -->

```
├── Figures_Individual
│   ├── 1_Preprocessing
│   │   ├── ATAC__peaks__average_profile
│   │   ├── ATAC__peaks__coverage
│   │   ├── ATAC__peaks__grouped_plots
│   │   ├── ATAC__peaks__saturation_curve
│   │   └── ATAC__reads__insert_size
│   ├── 2_Differential_Abundance
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
│   └── 3_Enrichment
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
│   ├── 2_Differential_Abundance
│   └── 3_Enrichment
├── Processed_Data
│   ├── 1_Preprocessing
│   │   ├── ATAC__peaks__annotated_rds
│   │   ├── ATAC__peaks__raw
│   │   ├── ATAC__peaks__split
│   │   ├── ATAC__peaks__split__no_BL
│   │   ├── ATAC__peaks__split__no_BL_input_control
│   │   ├── ATAC__peaks__split__no_BL_input_RNAi
│   │   ├── ATAC__reads__bam_asBed_atacShift_extendedBed
│   │   ├── ATAC__reads__bam_no_lowQ_dupli_mito
│   │   ├── ATAC__reads__fastqc_raw
│   │   ├── ATAC__reads__fastqc_trimmed
│   │   ├── ATAC__reads__multiQC
│   │   ├── mRNA__fastqc
│   │   ├── mRNA__kallisto_output
│   │   └── mRNA__multiQC
│   ├── 2_Differential_Abundance
│   │   ├── ATAC__all_peaks__bed
│   │   ├── ATAC__all_peaks__ChIPseeker
│   │   ├── ATAC__all_peaks__dataframe
│   │   ├── ATAC__all_peaks__DiffBind
│   │   ├── ATAC__all_peaks__gRange
│   │   ├── DA_split__bed_regions
│   │   ├── DA_split__genes_rds
│   │   ├── mRNA__all_genes__bed_promoters
│   │   ├── mRNA__all_genes__dataframe
│   │   └── mRNA__all_genes__rsleuth
│   ├── 3_Enrichment
│   │   └── motifs
│   └── ATAC__reads__bigwig_raw
├── Tables_Individual
│   ├── 1_Preprocessing
│   ├── 2_Differential_Abundance
│   │   └── ATAC_detailed
│   └── 3_Enrichment
│       ├── CHIP
│       ├── chrom_states
│       ├── func_anno_BP
│       ├── func_anno_KEGG
│       ├── genes_self
│       ├── motifs
│       └── peaks_self
└── Tables_Merged
    ├── 1_Preprocessing
    ├── 2_Differential_Abundance
    └── 3_Enrichment
```
