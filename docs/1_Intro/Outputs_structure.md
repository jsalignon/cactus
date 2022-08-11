

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [Data](/docs/2_Install/Data.md), [Test_datasets](/docs/2_Install/Test_datasets.md)
* [Run](/docs/3_Run/3_Run.md): [Input Data](/docs/3_Run/Input_data.md), [Input Files](/docs/3_Run/Input_files.md), [Parameters](/docs/3_Run/Parameters.md)

[](END_OF_MENU)


# Types

[Note: There should be one example Figure for each output]
There are three main types of outputs: processed data, figures and tables. Both figures and tables are available either merged (for easier scrolling) and separated (for lower memory usage). Tables can be either csv files or formatted Excel files.
Here is the structure of the output file with a short description of each output file and examples when relevant: 
TO DO. [add the tree of output as clickable and with a submenu, and description and image for each item]


# Structure

<!-- tree -I Run_Info -d -L 3 /home/jersal/workspace/cactus/test_datasets/human/results/25.07.22/ -->

```
.
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
