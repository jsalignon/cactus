
<img src="/docs/images/logo_cactus.png" width="400" />

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [References](/docs/2_Install/References.md), [Test datasets](/docs/2_Install/Test_datasets.md)
* [Inputs](/docs/3_Inputs/3_Inputs.md): [Data](/docs/3_Inputs/Data.md), [Design](/docs/3_Inputs/Design.md), [Parameters](/docs/3_Inputs/Parameters.md)
* [1. Preprocessing](/docs/4_Prepro/4_Prepro.md): [ATAC reads](/docs/4_Prepro/ATAC_reads.md), [ATAC peaks](/docs/4_Prepro/ATAC_peaks.md), [mRNA](/docs/4_Prepro/mRNA.md)
* [2. Differential Abundance](/docs/5_DA/5_DA.md): [ATAC](/docs/5_DA/DA_ATAC.md), [mRNA](/docs/5_DA/DA_mRNA.md), [Split](/docs/5_DA/Split.md)
* [3. Enrichment](/docs/6_Enrich/6_Enrich.md): [Enrichment](/docs/6_Enrich/Enrichment.md), [Figures](/docs/6_Enrich/Figures.md), [Tables](/docs/6_Enrich/Tables.md)

[](END_OF_MENU)



## Introduction

Tests datasets have beed generated for each of the 4 species supported by Cactus. This allows to make sure the pipeline works equally well for each species. Downloading the smallest dataset (fly) should be sufficient for users to get a sense on how the pipeline works. Here are the sizes of the test datasets:

| specie | raw_files | sampled_files | sampled_atac | sampled_mrna |  tar   |
|:------:|:---------:|:-------------:|:------------:|:------------:|:------:|
|  fly   |   17 GB   |    377 MB     |    329 MB    |    48 MB     | 360 MB |
|  worm  |   68 GB   |    1.4 GB     |    1.3 GB    |    41 MB     | 1.3 GB |
| mouse  |   52 GB   |     6 GB      |    5.9 GB    |    142 MB    | 5.7 GB |
| human  |  118 GB   |     9 GB      |    8.8 GB    |    190 MB    | 8.5 GB |

> **_NOTE:_**  Sampled mRNA-Seq datasets are similar between the 4 species, however, sampled ATAC-Seq datasets are much larger for human and mice than for worm and fly. This is due to large difference in genome sizes but not in transcriptome size between these species; as shown here:

| species |  genome | transcriptome | 
|:------:|:-------:|:-------------:|
|  fly   |  100 MB |     53 MB     |
|  worm  |  144 MB |     89 MB     |
| human  | 2731 MB |    158 MB     |
| mouse  | 3100 MB |    261 MB     |   

<!--  -->


<!-- 

get_markdown_table <- function(specie){
  tmp = read.table(paste0('~/workspace/cactus/test_datasets/preprocessing/', specie, '/samplesheet/samples_info_1.tsv'), header = T)
  knitr::kable(tmp, 'pipe', align = 'c')
}

get_markdown_table('worm')
get_markdown_table('fly')
get_markdown_table('human')
get_markdown_table('mouse')

setwd('~/workspace/cactus/test_datasets')

dir_size <- function(path, recursive = TRUE) {
  stopifnot(is.character(path))
  files <- list.files(path, full.names = T, recursive = recursive)
  vect_size <- sapply(files, function(x) file.size(x))
  size_files <- sum(vect_size)
  order = ifelse(size_files > 10^9, 'GB', 'MB')
  if(order == 'MB') size_files = paste0(format(size_files/10^6, digits = 2), ' MB')
  if(order == 'GB') size_files = paste0(format(size_files/10^9, digits = 2), ' GB')
  size_files
}

get_test_dataset_size <- function(specie){
  df = data.frame(
    specie = specie,
    raw_files = dir_size(paste0('preprocessing/', specie, '/fastq')),
    sampled_files =dir_size(paste0(specie, '/data')),
    sampled_atac = dir_size(paste0(specie, '/data/atac')),
    sampled_mrna = dir_size(paste0(specie, '/data/mrna')),
    stringsAsFactors = F
  )
  return(df)
}

df = do.call(rbind, lapply(c('fly', 'worm', 'mouse', 'human'), get_test_dataset_size))

df$tar = c('360 MB', '1.3 GB', '5.7 GB', '8.5 GB') 
# => see bash code below for the size of the tar files (to enter manually)

knitr::kable(df, 'pipe', align = 'c')

 -->


<!-- 

homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/test_datasets
refs_dir=$cactus_dir/references
arch_dir=$cactus_dir/figshare/files_to_upload
du -h $test_dir/**/data
ls -sh $arch_dir/*test.tar.gz

 -->

## Downloading test datasets

The test datasets can be downloaded with this command: 
```
nextflow run jsalignon/cactus/scripts/download/download.nf --test_datasets --species worm -r main -latest
```

The parameters for this command are:
 - *--species*: can be any of the 4 species supported by Cactus (worm, fly, mouse or human)
 - *--threads* can be set to determine the number of thread used by pigz for uncompressing the references archive files

 > **_NOTE:_**  The test datasets contains 4 folders: fastq, tsv, conf and yml; as described in the [Inputs](/docs/3_Inputs/3_Inputs.md) section.  
 
 > **_NOTE:_**  The test datasets should not be downloaded in the main cactus directory, otherwise the test dataset's conf folder will erase cactus' conf folder.  


##  Details on the test datasets origins and labels

<!-- 

homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
test_dir=$cactus_dir/test_datasets
cd $test_dir

R

library(magrittr)

get_kable_test_ds <- function(species){
  df = read.table(paste0('preprocessing/', species, '/samplesheet/samples_info_1.tsv'), header = T)
  colnames(df) = c('srx_id', 'srr_id', 'gsm_id', 'library_layout', 'library_strategy', 'sample_title', 'sample_id')
  df1 = df %>% dplyr::select(sample_id, library_strategy, sample_title, srr_id:library_layout)
  knitr::kable(df1, align = 'c')
}

get_kable_test_ds('worm')
get_kable_test_ds('human')
get_kable_test_ds('fly')
get_kable_test_ds('mouse')


 -->

## Worm and human ([GSE98758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98758))

  **GEO Status**: Public on Aug 29, 2018

  **GEO Title**: Genome-wide DNA accessibility maps and differential gene expression using ChIP-seq, ATAC-seq and RNA-seq for the human secondary fibroblast cell line hiF-T and whole worms with and without knockdown of FACT complex
  
  **GEO Summary**: To assess the mechanisms by which FACT depletion leads to increased sensitivity of cells to be reprogrammed, we measured the chromatin accessibility landscape using ATAC-seq following mock treatment, SSRP1 knockdown, or SUPT16H knockdown in human fibroblasts and mock, hmg-3 or hmg-4 knockdown in whole worms, and differential gene expression in hmg-3 knockout mutants or following mock, hmg-4, or spt-16 knockdown by RNAseq.
  
  **GEO Design**: Examination of two FACT complex components in human cells and worms with ChIP-seq, ATAC-seq and RNA-seq
  
  **Citation**: [*Kolundzic E, Ofenbauer A, Bulut SI, Uyar B et al. FACT Sets a Barrier for Cell Fate Reprogramming in Caenorhabditis elegans and Human Cells. Dev Cell 2018 Sep 10;46(5):611-626.e12. PMID: 30078731*](https://www.sciencedirect.com/science/article/pii/S1534580718305598)

  **Abstract**: The chromatin regulator FACT (facilitates chromatin transcription) is essential for ensuring stable gene expression by promoting transcription. In a genetic screen using Caenorhabditis elegans, we identified that FACT maintains cell identities and acts as a barrier for transcription factor-mediated cell fate reprogramming. Strikingly, FACT’s role as a barrier to cell fate conversion is conserved in humans as we show that FACT depletion enhances reprogramming of fibroblasts. Such activity is unexpected because FACT is known as a positive regulator of gene expression, and previously described reprogramming barriers typically repress gene expression. While FACT depletion in human fibroblasts results in decreased expression of many genes, a number of FACT-occupied genes, including reprogramming-promoting factors, show increased expression upon FACT depletion, suggesting a repressive function of FACT. Our findings identify FACT as a cellular reprogramming barrier in C. elegans and humans, revealing an evolutionarily conserved mechanism for cell fate protection.
  
  **Homologs**: 
|        Worm | Humans      |
|   :----:    |   :----:    |
|hmg-4, hmg-3 | SSRP1       |
|      spt-16 | SUPT16H     |

  **Worm samples**:

| sample_id | library_strategy | sample_title  |   srr_id   |   gsm_id   | library_layout |
|:---------:|:----------------:|:-------------:|:----------:|:----------:|:--------------:|
|   input   |     ATAC-seq     | input_control | SRR5000684 | GSM2385318 |     PAIRED     |
|   ctl_1   |     RNA-Seq      |   ctrl_rep1   | SRR5860412 | GSM2715402 |     SINGLE     |
|   ctl_2   |     RNA-Seq      |   ctrl_rep2   | SRR5860413 | GSM2715403 |     SINGLE     |
|   ctl_3   |     RNA-Seq      |   ctrl_rep3   | SRR5860414 | GSM2715404 |     SINGLE     |
|  hmg4_1   |     RNA-Seq      |   hmg4_rep1   | SRR5860415 | GSM2715405 |     SINGLE     |
|  hmg4_2   |     RNA-Seq      |   hmg4_rep2   | SRR5860416 | GSM2715406 |     SINGLE     |
|  hmg4_3   |     RNA-Seq      |   hmg4_rep3   | SRR5860417 | GSM2715407 |     SINGLE     |
|  spt16_1  |     RNA-Seq      |  spt16_rep1   | SRR5860418 | GSM2715408 |     SINGLE     |
|  spt16_2  |     RNA-Seq      |  spt16_rep2   | SRR5860419 | GSM2715409 |     SINGLE     |
|  spt16_3  |     RNA-Seq      |  spt16_rep3   | SRR5860420 | GSM2715410 |     SINGLE     |
|   ctl_1   |     ATAC-seq     |   Rluc_rep1   | SRR5860424 | GSM2715414 |     PAIRED     |
|   ctl_2   |     ATAC-seq     |   Rluc_rep2   | SRR5860425 | GSM2715415 |     PAIRED     |
|   ctl_3   |     ATAC-seq     |   Rluc_rep3   | SRR5860426 | GSM2715416 |     PAIRED     |
|  spt16_1  |     ATAC-seq     |  spt-16_rep1  | SRR5860430 | GSM2715420 |     PAIRED     |
|  spt16_2  |     ATAC-seq     |  spt-16_rep2  | SRR5860431 | GSM2715421 |     PAIRED     |
|  spt16_3  |     ATAC-seq     |  spt-16_rep3  | SRR5860432 | GSM2715422 |     PAIRED     |
|  hmg4_1   |     ATAC-seq     |  hmg-4_rep1   | SRR5860433 | GSM2715423 |     PAIRED     |
|  hmg4_2   |     ATAC-seq     |  hmg-4_rep2   | SRR5860434 | GSM2715424 |     PAIRED     |
|  hmg4_3   |     ATAC-seq     |  hmg-4_rep3   | SRR5860435 | GSM2715425 |     PAIRED     |


  **Human samples**:
  
| sample_id | library_strategy |                   sample_title                    |   srr_id   |   gsm_id   | library_layout |
|:---------:|:----------------:|:-------------------------------------------------:|:----------:|:----------:|:--------------:|
|   ctl_1   |     ATAC-seq     |                mock_ATAC-seq_rep1                 | SRR5521292 | GSM2611319 |     PAIRED     |
|   ctl_2   |     ATAC-seq     |                mock_ATAC-seq_rep2                 | SRR5521293 | GSM2611320 |     PAIRED     |
|  ssrp1_1  |     ATAC-seq     |                SSRP1_ATAC-seq_rep1                | SRR5521294 | GSM2611321 |     PAIRED     |
|  ssrp1_2  |     ATAC-seq     |                SSRP1_ATAC-seq_rep2                | SRR5521295 | GSM2611322 |     PAIRED     |
| supt16h_1 |     ATAC-seq     |               SUPT16H_ATAC-seq_rep1               | SRR5521296 | GSM2611323 |     PAIRED     |
| supt16h_2 |     ATAC-seq     |               SUPT16H_ATAC-seq_rep2               | SRR5521297 | GSM2611324 |     PAIRED     |
|   ctl_1   |     RNA-Seq      |     pri_mockTotal_A:_human_mock_RNA-seq_rep1      | SRR7101006 | GSM3127942 |     PAIRED     |
|   ctl_2   |     RNA-Seq      |     pri_mockTotal_B:_human_mock_RNA-seq_rep2      | SRR7101007 | GSM3127943 |     PAIRED     |
|   ctl_3   |     RNA-Seq      |     pri_mockTotal_C:_human_mock_RNA-seq_rep3      | SRR7101008 | GSM3127944 |     PAIRED     |
|  ssrp1_1  |     RNA-Seq      |   pri_ssrp1Total_A:_human_Ssrp1_kd_RNA-seq_rep1   | SRR7101009 | GSM3127945 |     PAIRED     |
|  ssrp1_2  |     RNA-Seq      |   pri_ssrp1Total_B:_human_Ssrp1_kd_RNA-seq_rep2   | SRR7101010 | GSM3127946 |     PAIRED     |
|  ssrp1_3  |     RNA-Seq      |   pri_ssrp1Total_C:_human_Ssrp1_kd_RNA-seq_rep3   | SRR7101011 | GSM3127947 |     PAIRED     |
| supt16h_1 |     RNA-Seq      | pri_supt16hTotal_A:_human_Supt16h_kd_RNA-seq_rep1 | SRR7101012 | GSM3127948 |     PAIRED     |
| supt16h_2 |     RNA-Seq      | pri_supt16hTotal_B:_human_Supt16h_kd_RNA-seq_rep2 | SRR7101013 | GSM3127949 |     PAIRED     |
| supt16h_3 |     RNA-Seq      | pri_supt16hTotal_C:_human_Supt16h_kd_RNA-seq_rep3 | SRR7101014 | GSM3127950 |     PAIRED     |


## Fly ([GSE149339](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149336))


  **GEO Status**: Public on May 10, 2020

  **GEO Title**: Pioneer factor GAF cooperates with PBAP and NURF to regulate transcription
  
  **GEO Summary**: The Drosophila pioneer factor GAF is known to be essential for RNA Pol II promoter-proximal pausing and the removal of nucleosomes from a set of target promoters with GAGAG motifs. We and others have speculated that GAF recruits the ISWI family ATP-dependent chromatin remodeling complex NURF, on the basis that NURF and GAF are both required to remodel nucleosomes on an hsp70 promoter in vitro and that GAF interacts physically with NURF. However, GAF was also recently shown to interact with PBAP, a SWI/SNF family remodeler. To test which of these remodeling complexes GAF works with, we depleted GAF, NURF301, BAP170, and NURF301+BAP170 in Drosophila S2 cells using RNAi. We used a combination of PRO-seq, ATAC-seq, 3'RNA-seq, and CUT&RUN to demonstrate that while GAF and PBAP synergistically open chromatin at target promoters which allows Pol II recruitment and pausing to proceed, GAF and NURF also synergistically position the +1 nucleosome to ensure efficient pause release and transition to productive elongation.
  
  **GEO Design**: We treated two independent replicates of Drosophila S2 cells with dsRNA to LACZ (control), GAF, NURF301, BAP170 (the unique subunits of the NURF and PBAP complexes, respectively), and NURF301+BAP170. After 5 days, we harvested cells, validated knockdowns, and performed PRO-seq, ATAC-seq and 3'RNA-seq. We also performed CUT&RUN for both GAF and NURF301 in untreated S2 cells.
  
  **Citation**:  [*Judd, J., Duarte, F. M. & Lis, J. T. Pioneer-like factor GAF cooperates with PBAP (SWI/SNF) and NURF (ISWI) to regulate transcription. Genes Dev. 35, 147–156 (2021).*](http://genesdev.cshlp.org/content/35/1-2/147)

  **Abstract**: Transcriptionally silent genes must be activated throughout development. This requires nucleosomes be removed from promoters and enhancers to allow transcription factor (TF) binding and recruitment of coactivators and RNA polymerase II (Pol II). Specialized pioneer TFs bind nucleosome-wrapped DNA to perform this chromatin opening by mechanisms that remain incompletely understood. Here, we show that GAGA factor (GAF), a Drosophila pioneer-like factor, functions with both SWI/SNF and ISWI family chromatin remodelers to allow recruitment of Pol II and entry to a promoter-proximal paused state, and also to promote Pol II's transition to productive elongation. We found that GAF interacts with PBAP (SWI/SNF) to open chromatin and allow Pol II to be recruited. Importantly, this activity is not dependent on NURF as previously proposed; however, GAF also synergizes with NURF downstream from this process to ensure efficient Pol II pause release and transition to productive elongation, apparently through its role in precisely positioning the +1 nucleosome. These results demonstrate how a single sequence-specific pioneer TF can synergize with remodelers to activate sets of genes. Furthermore, this behavior of remodelers is consistent with findings in yeast and mice, and likely represents general, conserved mechanisms found throughout eukarya. 
  
  **Samples**:
  
| sample_id  | library_strategy |        sample_title        |   srr_id    |   gsm_id   | library_layout |
|:----------:|:----------------:|:--------------------------:|:-----------:|:----------:|:--------------:|
|   ctl_1    |     ATAC-seq     |     LACZ_ATACseq_Rep1      | SRR11607688 | GSM4498282 |     PAIRED     |
|   ctl_1    |     ATAC-seq     |     LACZ_ATACseq_Rep1      | SRR11607689 | GSM4498282 |     PAIRED     |
|   ctl_2    |     ATAC-seq     |     LACZ_ATACseq_Rep2      | SRR11607690 | GSM4498283 |     PAIRED     |
|   ctl_2    |     ATAC-seq     |     LACZ_ATACseq_Rep2      | SRR11607691 | GSM4498283 |     PAIRED     |
|   gaf_1    |     ATAC-seq     |      GAF_ATACseq_Rep1      | SRR11607692 | GSM4498284 |     PAIRED     |
|   gaf_1    |     ATAC-seq     |      GAF_ATACseq_Rep1      | SRR11607693 | GSM4498284 |     PAIRED     |
|   gaf_2    |     ATAC-seq     |      GAF_ATACseq_Rep2      | SRR11607675 | GSM4498285 |     PAIRED     |
|   gaf_2    |     ATAC-seq     |      GAF_ATACseq_Rep2      | SRR11607694 | GSM4498285 |     PAIRED     |
|   b170_1   |     ATAC-seq     |    BAP170_ATACseq_Rep1     | SRR11607676 | GSM4498286 |     PAIRED     |
|   b170_1   |     ATAC-seq     |    BAP170_ATACseq_Rep1     | SRR11607677 | GSM4498286 |     PAIRED     |
|   b170_2   |     ATAC-seq     |    BAP170_ATACseq_Rep2     | SRR11607678 | GSM4498287 |     PAIRED     |
|   b170_2   |     ATAC-seq     |    BAP170_ATACseq_Rep2     | SRR11607679 | GSM4498287 |     PAIRED     |
|   n301_1   |     ATAC-seq     |    NURF301_ATACseq_Rep1    | SRR11607680 | GSM4498288 |     PAIRED     |
|   n301_1   |     ATAC-seq     |    NURF301_ATACseq_Rep1    | SRR11607681 | GSM4498288 |     PAIRED     |
|   n301_2   |     ATAC-seq     |    NURF301_ATACseq_Rep2    | SRR11607682 | GSM4498289 |     PAIRED     |
|   n301_2   |     ATAC-seq     |    NURF301_ATACseq_Rep2    | SRR11607683 | GSM4498289 |     PAIRED     |
| n301b170_1 |     ATAC-seq     | NURF301BAP170_ATACseq_Rep1 | SRR11607684 | GSM4498290 |     PAIRED     |
| n301b170_1 |     ATAC-seq     | NURF301BAP170_ATACseq_Rep1 | SRR11607685 | GSM4498290 |     PAIRED     |
| n301b170_2 |     ATAC-seq     | NURF301BAP170_ATACseq_Rep2 | SRR11607686 | GSM4498291 |     PAIRED     |
| n301b170_2 |     ATAC-seq     | NURF301BAP170_ATACseq_Rep2 | SRR11607687 | GSM4498291 |     PAIRED     |
|   gaf_2    |     RNA-Seq      |      GAF_RNAseq_Rep2       | SRR11607698 | GSM4498295 |     SINGLE     |
|   b170_1   |     RNA-Seq      |     BAP170_RNAseq_Rep1     | SRR11607699 | GSM4498296 |     SINGLE     |
|   b170_2   |     RNA-Seq      |     BAP170_RNAseq_Rep2     | SRR11607700 | GSM4498297 |     SINGLE     |
|   n301_1   |     RNA-Seq      |    NURF301_RNAseq_Rep1     | SRR11607701 | GSM4498298 |     SINGLE     |
|   n301_2   |     RNA-Seq      |    NURF301_RNAseq_Rep2     | SRR11607702 | GSM4498299 |     SINGLE     |
| n301b170_1 |     RNA-Seq      | NURF301BAP170_RNAseq_Rep1  | SRR11607703 | GSM4498300 |     SINGLE     |
| n301b170_2 |     RNA-Seq      | NURF301BAP170_RNAseq_Rep2  | SRR11607704 | GSM4498301 |     SINGLE     |
|   ctl_1    |     RNA-Seq      |      LACZ_RNAseq_Rep1      | SRR11607695 | GSM4498292 |     SINGLE     |
|   ctl_2    |     RNA-Seq      |      LACZ_RNAseq_Rep2      | SRR11607696 | GSM4498293 |     SINGLE     |
|   gaf_1    |     RNA-Seq      |      GAF_RNAseq_Rep1       | SRR11607697 | GSM4498294 |     SINGLE     |


## Mouse ([GSE193393](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193393))


  **GEO Status**: Public on Jun 23, 2022

  **GEO Title**: PHF20 Activates Autophagy Genes through Enhancer Activation via H3K36me2 Binding Activity
  
  **GEO Summary**: Autophagy is a catabolic pathway that maintains cellular homeostasis under various stress conditions, including nutrient-deprived conditions. To elevate autophagic flux to a sufficient level under stress conditions, transcriptional activation of autophagy genes occurs to replenish autophagy components. Here, using combination of RNA-seq, ATAC-seq and ChIP-seq, we demonstrated found that plant homeodomain finger protein 20 (Phf20PHF20), which is an epigenetic reader possessing methyl binding activity, plays a key role in controlling the expression of autophagy genes. PHF20 activates autophagy genes through enhancer activation via H3K36me2 binding activity as an epigenetic reader and that our findings emphasize the importance of nuclear regulation of autophagy.
  
  **GEO Design**: mRNA-seq, ATAC-seq and ChIP-seq experiments under normal and 24hrs of glucose starvation condition in WT and Phf20-/- MEFs
  
  **Citation**: Park SW, Kim J, Oh S, Lee J et al. PHF20 is crucial for epigenetic control of starvation-induced autophagy through enhancer activation. Nucleic Acids Res 2022 Aug 12;50(14):7856-7872. PMID: 35821310
  
  **Samples**:
  
| sample_id  | library_strategy |          sample_title           |   srr_id    |   gsm_id   | library_layout |
|:----------:|:----------------:|:-------------------------------:|:-----------:|:----------:|:--------------:|
|    Wt_1    |     ATAC-seq     |    WT_control_ATAC-seq_rep1     | SRR17483668 | GSM5776726 |     PAIRED     |
|    Wt_2    |     ATAC-seq     |    WT_control_ATAC-seq_rep2     | SRR17483667 | GSM5776727 |     PAIRED     |
| WtStarv_1  |     ATAC-seq     |    WT_GlcStarv_ATAC-seq_rep1    | SRR17483666 | GSM5776728 |     PAIRED     |
| WtStarv_2  |     ATAC-seq     |    WT_GlcStarv_ATAC-seq_rep2    | SRR17483665 | GSM5776729 |     PAIRED     |
|   Phf_1    |     ATAC-seq     | Phf20-/-_control_ATAC-seq_rep1  | SRR17483664 | GSM5776730 |     PAIRED     |
|   Phf_2    |     ATAC-seq     | Phf20-/-_control_ATAC-seq_rep2  | SRR17483663 | GSM5776731 |     PAIRED     |
| PhfStarv_1 |     ATAC-seq     | Phf20-/-_GlcStarv_ATAC-seq_rep1 | SRR17483662 | GSM5776732 |     PAIRED     |
| PhfStarv_2 |     ATAC-seq     | Phf20-/-_GlcStarv_ATAC-seq_rep2 | SRR17483661 | GSM5776733 |     PAIRED     |
|    Wt_1    |     RNA-Seq      |     WT_control_RNA-seq_rep1     | SRR17535397 | GSM5799507 |     PAIRED     |
|    Wt_2    |     RNA-Seq      |     WT_control_RNA-seq_rep2     | SRR17535396 | GSM5799508 |     PAIRED     |
| WtStarv_1  |     RNA-Seq      |    WT_GlcStarv_RNA-seq_rep1     | SRR17535395 | GSM5799509 |     PAIRED     |
| WtStarv_2  |     RNA-Seq      |    WT_GlcStarv_RNA-seq_rep2     | SRR17535394 | GSM5799510 |     PAIRED     |
|   Phf_1    |     RNA-Seq      |  Phf20-/-_control_RNA-seq_rep1  | SRR17535392 | GSM5799511 |     PAIRED     |
|   Phf_2    |     RNA-Seq      |  Phf20-/-_control_RNA-seq_rep2  | SRR17535391 | GSM5799512 |     PAIRED     |
| PhfStarv_1 |     RNA-Seq      | Phf20-/-_GlcStarv_RNA-seq_rep1  | SRR17535390 | GSM5799513 |     PAIRED     |
| PhfStarv_2 |     RNA-Seq      | Phf20-/-_GlcStarv_RNA-seq_rep2  | SRR17535393 | GSM5799514 |     PAIRED     |

