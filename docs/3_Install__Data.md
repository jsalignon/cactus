
* Introduction: [Cactus](/README.md), [Quick Start](/docs/1_Intro__Quick_start.md), 
* Overview: [Graph](/docs/2_Overview__Graph.md), [Outputs structure](/docs/2_Overview__Outputs_structure.md)
* Install: [Dependencies](/docs/3_Install__Dependencies.md), [Containers](/docs/3_Install__Containers.md), [Data](/docs/3_Install__Data.md), [Test_datasets](/docs/3_Install__Test_datasets.md)
* Run: [Input Data](/docs/4_Run__Input_data.md), [Input Files](/docs/4_Run__Input_files.md), [Parameters](/docs/4_Run__Parameters.md)
* Preprocessing: ATAC: [Method](/docs/5_AP__Method.md), [Figures](/docs/5_AP__Figures.md), [MultiQC](/docs/5_AP__MultiQC.md), mRNA: [Method](/docs/6_MP__Method.md), [MultiQC](/docs/6_MP__MultiQC.md)
* Differential Abundance: [ATAC](/docs/7_DA__DiffBind.md), [mRNA](/docs/7_DA__Sleuth.md), [Figures](/docs/7_DA__Figures.md), [Tables](/docs/7_DA__Tables.md)
* Splitting peak sets: [Split](/docs/8_SP__Splitting_peak_sets.md), [Venn Diagrams](/docs/8_SP__Venn_diagrams.md)
* Enrichment: [ChIP](/docs/9_Enrich__CHIP.md), [Chromatin state](/docs/9_Enrich__Chromatin_states.md), [Motifs](/docs/9_Enrich__Motifs.md), [Func. Anno.](/docs/9_Enrich__Functional_annotations.md), [Figures](/docs/9_Enrich__Figures.md), [Tables](/docs/9_Enrich__Tables.md)


[](END_OF_MENU)

## Introduction
<!-- 
# get data folders sizes: du -h -d1 human mouse fly worm
# print table (in R):

tb = tibble::tribble(
  ~Specie, ~Size,  ~Assembly, ~Nickname, ~Ensembl_Release,
   'worm',   1.3, 'WBcel235',    'ce11',    '107',
    'fly',   1.5,    'BDGP6',     'dm6',    '107',
  'mouse',    12,   'GRCm38',    'mm10',    '102',
  'human',    20,   'GRCh38',    'hg38',    '107'
)
knitr::kable(tb, 'pipe', align = c('c', 'r', 'c', 'c', 'c'))

-->

Cactus uses a set of preparsed references that are used throughout the scripts. These are automatically downloaded the first time the pipeline is run for the given specie under investigation. The size of the Cactus data for each specie is:

| Specie | Size| Assembly | Nickname | Ensembl_Release |
|:------:|----:|:--------:|:--------:|:---------------:|
|  worm  |  1.3| WBcel235 |   ce11   |       107       |
|  fly   |  1.5|  BDGP6   |   dm6    |       107       |
| mouse  | 12.0|  GRCm38  |   mm10   |       102       |
| human  | 20.0|  GRCh38  |   hg38   |       107       |

The latest Ensembl release has been used for worm, fly and mouse (July 2022). For mouse, the older Ensembl release 102 (Nov. 2020) has been used together with the genome assembly mm10, since mm39 is not yet available in Homer and in the Encode bed files.

Note that in order to save space, on can download the specie dataset, and then keep only the chromatin states of interests. Indeed, only one chromatin state file is used by Cactus, but several are available for users to chose depending on their conditions under study. For mouse and humans this should reduce the size of the dataset by respectively about 8 Gb and 1 Gb.


## Pipeline to get data

![DAG](/docs/images/dag_get_data.svg)

The pipeline download the data from various sources and process it to produce all the files required by cactus. These can be grouped in 6 broad categories: parsed genome sequences, parsed genome annotations, CHIP-Seq files, chromatin states files, parsed motifs files and software (Homer data), bowtie2 indexes of a contaminant genome. In order to maxime reproducibility, all processes (excepting the one to download HiHMM files) are encapsulated within containers. Here are some details on what the parsing pipeline does:
  - Motifs:
   - [Homer data](http://homer.ucsd.edu/homer/): The latest genome (6.4), organism (6.3) and promoters (5.5) versions were downloaded from http://homer.ucsd.edu/homer/data
   - [CISBP Motifs](http://cisbp.ccbr.utoronto.ca/): The TF_Information.txt file is downloaded from http://cisbp.ccbr.utoronto.ca/data. This file contains for each transcription factor, either direct motifs, or if not available, motifs inferred with more than 90% confidence (all motifs with the best score are kept). There are 10,329 motifs from the Encod species. 2,779 motifs were kept after removing duplicated motifs (inferred motifs with the same score), motifs with empty names and empty motifs (worm: 372, fly: 424, mouse: 872, human: 1,111). Homer threshold were computed as: $0.5 x motif_log2_odd_score$, with $motif_log2_odd_score = sum(log2(max_weight_by_position / 0.25))$. Finally, for each specie, a homer_motifs.txt file is generated that contains consensus sequence, TF name, homer threshold, and the motif itself.
   
  - genome sequence: 





## Structure of the data folders
<!-- tree -d -L 3 --filelimit 5 human -->

```
.
├── bowtie2_indexes_conta
├── CHIP
├── chromatin_states
├── genome
│   ├── annotation
│   │   ├── bed_regions
│   │   ├── filtered
│   │   └── R
│   └── sequence
│       └── bowtie2_indexes
└── homer_data
    ├── accession
    ├── genomes
    │   └── /nickname/
    ├── GO
    └── promoters
```
