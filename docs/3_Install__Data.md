
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
    - [Homer data](https://doi.org/10.1016/j.molcel.2010.05.004): The latest genome (6.4), organism (6.3) and promoters (5.5) versions were downloaded from http://homer.ucsd.edu/homer.
    - [CISBP Motifs](https://doi.org/10.1016/j.cell.2014.08.009): The TF_Information.txt file is downloaded from http://cisbp.ccbr.utoronto.ca. This file contains for each transcription factor, either direct motifs, or if not available, motifs inferred with more than 90% confidence (all motifs with the best score are kept). There are 10,329 motifs from the Encod species. 2,779 motifs were kept after removing duplicated motifs (inferred motifs with the same score), motifs with empty names and empty motifs (worm: 372, fly: 424, mouse: 872, human: 1,111). Homer threshold were computed as: *0.5 x motif_log2_odd_score*, with *motif_log2_odd_score = sum(log2(max_weight_by_position / 0.25))*. Finally, for each specie, a homer_motifs.txt file is generated that contains consensus sequence, TF name, homer threshold, and the motif itself.
   
  - [Blacklisted regions](https://doi.org/10.1038/s41598-019-45839-z): Blacklisted regions were donwloaded from https://github.com/Boyle-Lab/Blacklist/tree/master/lists. Contig names were shifted from NCBI to Ensembl names using the UpdateContigNames function from [cvbio](https://github.com/clintval/cvbio).
  
  <!--  head -2 **/available_chip_ontology_groups.txt -->
  - CHIP-Seq: The [ENCODE API](https://www.encodeproject.org/help/rest-api/) was used to get data and metadata. 2,714 CHIP-Seq bed files were selected and downloaded (worm: 473, fly: 531, mouse: 156, human: 1,554) using these filters: assay_title = "TF ChIP-seq" and output_type = "optimal IDR thresholded peaks". The slim annotations (cell, organ, development and system) were parsed and used to create groups of CHIP-Seq that share the same annotations and can be used for more detailled analysis (see [CHIP-Seq](/docs/3_Install__Data.md##CHIP-Seq) and [Parameters](/docs/4_Run__Parameters.md))
  -  profiles that can be used 
  -  




  
  1554
    - 2,716
  
  - genome: 
    - Sequence:
    - Annotations:


## CHIP-Seq


<!-- 
R
setwd('/home/jersal/workspace/cactus/data')

print_ontology_table <- function(specie){
  dt = data.table::fread(paste0(specie, '/available_chip_ontology_groups.txt'))[, 2:3]
  knitr::kable(dt, 'pipe', align = 'c', row.names = F)
}
v_species = c('worm', 'fly', 'mouse', 'human')
sapply(v_species, print_ontology_table)

-->

### Worm

#### Legend
-  L1: L1 stage
-  YA: Young Adult
-  L4: L4 stage
-  LE: Late Embryo
-  L3: L3 stage
-  L2: L2 stage
-  MS: Mixed Stage embryo
-  EE: Early Embryo
-  ME: Mid-Embryo
-  Da: Dauer

#### Distribution
```
	L1 YA L4 LE L3 L2 MS EE ME Da
	93 86 77 67 50 50 35 12  2  1
```

#### Ontology groups

|         ontology         | number_of_chip |
|:------------------------:|:--------------:|
| cell_line.whole_organism |      473       |
|           all            |      473       |


### Fly

#### Legend
- Em: Embryo
- PP: Prepupa
- WT: Wandering third instar larval stage
- Pu: Puppa
- FA: Female Adult
- WP: White prepupa stage
- MA: Male Adult
- Kc: Kc167 cell line (embryonic)
- SO: Cell line from strain Oregon-R S2
- MS: Mixed Sex Adult
- OF: Ovary female
- La: larva (48 days)

#### Distribution
```
	 Em  PP  WT  Pu  FA  WP  MA  Kc  SO  MS  OF  La
	372  57  46  14  10   8   8   7   3   3   1   1
```

#### Ontology groups

|         ontology         | number_of_chip |
|:------------------------:|:--------------:|
|           all            |      531       |
| cell_line.whole_organism |      520       |


### Mouse

#### Legend
- Other: Any ontology present less than 4 times
- CH12LX: CH12.LX
- G1EER4: G1E-ER4
-  ESE14: ES-E14

#### Distribution
```
   MEL CH12LX  Other  liver   lung  heart G1EER4    G1E  ESE14
    49     39     36      7      5      5      5      5      5
```

#### Ontology groups

|           ontology           | number_of_chip |
|:----------------------------:|:--------------:|
|             all              |      156       |
|     development.mesoderm     |      107       |
|     system.immune_system     |       96       |
|   system.digestive_system    |       62       |
|      organ.immune_organ      |       52       |
|         organ.spleen         |       50       |
|        cell_line.MEL         |       49       |
|    cell_type.cancer_cell     |       49       |
| cell_type.hematopoietic_cell |       43       |
|     cell_type.leukocyte      |       41       |
|      cell_line.CH12.LX       |       39       |
|       cell_type.B_cell       |       39       |
|       organ.lymph_node       |       39       |
|         organ.blood          |       39       |
|      organ.bodily_fluid      |       39       |
|         organ.embryo         |       23       |
|     cell_type.stem_cell      |       21       |
|     development.endoderm     |       19       |
|   cell_type.embryonic_cell   |       11       |


### Human

#### Legend
- NC: Neural Cell
- Other: Any ontology present less than 5 times
- Ishikaw: Ishikawa
- SKNSH: SK-N-SH

#### Distribution
```
   K562   HepG2  HEK293 GM12878   Other    MCF7    A549      H1   liver Ishikaw
    432     249     193     156     137     106      55      54      35      24

  SKNSH  HeLaS3 HEK293T   IMR90  MCF10A  HCT116    T47D GM12891      NC GM23338
     19      19      17      12      10      10       7       7       6       6
 ```

#### Ontology groups

|             ontology             | number_of_chip |
|:--------------------------------:|:--------------:|
|               all                |      1554      |
|      cell_type.cancer_cell       |      956       |
|       development.mesoderm       |      932       |
|    cell_type.epithelial_cell     |      645       |
|         organ.epithelium         |      645       |
|       system.immune_system       |      637       |
|           organ.blood            |      631       |
|   cell_type.hematopoietic_cell   |      630       |
|       cell_type.leukocyte        |      628       |
|        organ.bodily_fluid        |      624       |
|          cell_line.K562          |      432       |
|      system.exocrine_system      |      411       |
|       organ.exocrine_gland       |      411       |
|       development.endoderm       |      401       |
|     system.digestive_system      |      306       |
|     system.endocrine_system      |      289       |
|           organ.liver            |      284       |
|      organ.endocrine_gland       |      284       |
|         cell_line.HepG2          |      249       |
|     system.excretory_system      |      212       |
|           organ.kidney           |      212       |
|         cell_line.HEK293         |      210       |
|         cell_type.B_cell         |      188       |
|       development.ectoderm       |      187       |
|        cell_line.GM12878         |      156       |
|    system.integumental_system    |      143       |
|       organ.mammary_gland        |      127       |
|          cell_line.MCF7          |      106       |
|    system.respiratory_system     |       75       |
|            organ.lung            |       75       |
|    system.reproductive_system    |       71       |
|       cell_type.stem_cell        |       60       |
|           organ.embryo           |       59       |
|          cell_line.A549          |       55       |
|     cell_type.embryonic_cell     |       55       |
|           cell_line.H1           |       54       |
|  system.central_nervous_system   |       46       |
|           organ.uterus           |       43       |
|           organ.brain            |       38       |
|     organ.connective_tissue      |       36       |
|         cell_line.liver          |       35       |
|       cell_type.fibroblast       |       35       |
| cell_type.connective_tissue_cell |       32       |
|        cell_line.Ishikawa        |       24       |
|   cell_type.neuroblastoma_cell   |       20       |
|         cell_line.HeLaS3         |       19       |
|         cell_line.SKNSH          |       19       |
|      cell_type.neural_cell       |       19       |
|        cell_line.HEK293T         |       17       |
|        organ.skin_of_body        |       17       |
|      organ.large_intestine       |       15       |
|           organ.colon            |       15       |
|         organ.intestine          |       15       |
|       organ.prostate_gland       |       14       |
|         cell_line.IMR90          |       12       |
|    system.circulatory_system     |       11       |



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
