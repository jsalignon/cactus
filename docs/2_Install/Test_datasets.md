
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

get_kable_test_ds <- function(species){
  df = read.table(paste0('preprocessing/', species, '/samplesheet/samples_info_1.tsv'), header = T)
  knitr::kable(df[, c(1:5, 7)], align = 'c')
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

| experiment_accession | run_accession | sample_alias | library_layout | library_strategy | sample_title  | sample_id |
|:--------------------:|:-------------:|:------------:|:--------------:|:----------------:|:-------------:|:---------:|
|      SRX2333004      |  SRR5000684   |  GSM2385318  |     PAIRED     |     ATAC-seq     | input_control |   input   |
|      SRX3029112      |  SRR5860412   |  GSM2715402  |     SINGLE     |     RNA-Seq      |   ctrl_rep1   |   ctl_1   |
|      SRX3029113      |  SRR5860413   |  GSM2715403  |     SINGLE     |     RNA-Seq      |   ctrl_rep2   |   ctl_2   |
|      SRX3029114      |  SRR5860414   |  GSM2715404  |     SINGLE     |     RNA-Seq      |   ctrl_rep3   |   ctl_3   |
|      SRX3029115      |  SRR5860415   |  GSM2715405  |     SINGLE     |     RNA-Seq      |   hmg4_rep1   |  hmg4_1   |
|      SRX3029116      |  SRR5860416   |  GSM2715406  |     SINGLE     |     RNA-Seq      |   hmg4_rep2   |  hmg4_2   |
|      SRX3029117      |  SRR5860417   |  GSM2715407  |     SINGLE     |     RNA-Seq      |   hmg4_rep3   |  hmg4_3   |
|      SRX3029118      |  SRR5860418   |  GSM2715408  |     SINGLE     |     RNA-Seq      |  spt16_rep1   |  spt16_1  |
|      SRX3029119      |  SRR5860419   |  GSM2715409  |     SINGLE     |     RNA-Seq      |  spt16_rep2   |  spt16_2  |
|      SRX3029120      |  SRR5860420   |  GSM2715410  |     SINGLE     |     RNA-Seq      |  spt16_rep3   |  spt16_3  |
|      SRX3029124      |  SRR5860424   |  GSM2715414  |     PAIRED     |     ATAC-seq     |   Rluc_rep1   |   ctl_1   |
|      SRX3029125      |  SRR5860425   |  GSM2715415  |     PAIRED     |     ATAC-seq     |   Rluc_rep2   |   ctl_2   |
|      SRX3029126      |  SRR5860426   |  GSM2715416  |     PAIRED     |     ATAC-seq     |   Rluc_rep3   |   ctl_3   |
|      SRX3029130      |  SRR5860430   |  GSM2715420  |     PAIRED     |     ATAC-seq     |  spt-16_rep1  |  spt16_1  |
|      SRX3029131      |  SRR5860431   |  GSM2715421  |     PAIRED     |     ATAC-seq     |  spt-16_rep2  |  spt16_2  |
|      SRX3029132      |  SRR5860432   |  GSM2715422  |     PAIRED     |     ATAC-seq     |  spt-16_rep3  |  spt16_3  |
|      SRX3029133      |  SRR5860433   |  GSM2715423  |     PAIRED     |     ATAC-seq     |  hmg-4_rep1   |  hmg4_1   |
|      SRX3029134      |  SRR5860434   |  GSM2715424  |     PAIRED     |     ATAC-seq     |  hmg-4_rep2   |  hmg4_2   |
  |      SRX3029135      |  SRR5860435   |  GSM2715425  |     PAIRED     |     ATAC-seq     |  hmg-4_rep3   |  hmg4_3   |

  **Human samples**:
  
| experiment_accession | run_accession | sample_alias | library_layout | library_strategy |                   sample_title                    | sample_id |
|:--------------------:|:-------------:|:------------:|:--------------:|:----------------:|:-------------------------------------------------:|:---------:|
|      SRX2794533      |  SRR5521292   |  GSM2611319  |     PAIRED     |     ATAC-seq     |                mock_ATAC-seq_rep1                 |   ctl_1   |
|      SRX2794534      |  SRR5521293   |  GSM2611320  |     PAIRED     |     ATAC-seq     |                mock_ATAC-seq_rep2                 |   ctl_2   |
|      SRX2794535      |  SRR5521294   |  GSM2611321  |     PAIRED     |     ATAC-seq     |                SSRP1_ATAC-seq_rep1                |  ssrp1_1  |
|      SRX2794536      |  SRR5521295   |  GSM2611322  |     PAIRED     |     ATAC-seq     |                SSRP1_ATAC-seq_rep2                |  ssrp1_2  |
|      SRX2794537      |  SRR5521296   |  GSM2611323  |     PAIRED     |     ATAC-seq     |               SUPT16H_ATAC-seq_rep1               | supt16h_1 |
|      SRX2794538      |  SRR5521297   |  GSM2611324  |     PAIRED     |     ATAC-seq     |               SUPT16H_ATAC-seq_rep2               | supt16h_2 |
|      SRX4029346      |  SRR7101006   |  GSM3127942  |     PAIRED     |     RNA-Seq      |     pri_mockTotal_A:_human_mock_RNA-seq_rep1      |   ctl_1   |
|      SRX4029347      |  SRR7101007   |  GSM3127943  |     PAIRED     |     RNA-Seq      |     pri_mockTotal_B:_human_mock_RNA-seq_rep2      |   ctl_2   |
|      SRX4029348      |  SRR7101008   |  GSM3127944  |     PAIRED     |     RNA-Seq      |     pri_mockTotal_C:_human_mock_RNA-seq_rep3      |   ctl_3   |
|      SRX4029349      |  SRR7101009   |  GSM3127945  |     PAIRED     |     RNA-Seq      |   pri_ssrp1Total_A:_human_Ssrp1_kd_RNA-seq_rep1   |  ssrp1_1  |
|      SRX4029350      |  SRR7101010   |  GSM3127946  |     PAIRED     |     RNA-Seq      |   pri_ssrp1Total_B:_human_Ssrp1_kd_RNA-seq_rep2   |  ssrp1_2  |
|      SRX4029352      |  SRR7101011   |  GSM3127947  |     PAIRED     |     RNA-Seq      |   pri_ssrp1Total_C:_human_Ssrp1_kd_RNA-seq_rep3   |  ssrp1_3  |
|      SRX4029353      |  SRR7101012   |  GSM3127948  |     PAIRED     |     RNA-Seq      | pri_supt16hTotal_A:_human_Supt16h_kd_RNA-seq_rep1 | supt16h_1 |
|      SRX4029354      |  SRR7101013   |  GSM3127949  |     PAIRED     |     RNA-Seq      | pri_supt16hTotal_B:_human_Supt16h_kd_RNA-seq_rep2 | supt16h_2 |
|      SRX4029355      |  SRR7101014   |  GSM3127950  |     PAIRED     |     RNA-Seq      | pri_supt16hTotal_C:_human_Supt16h_kd_RNA-seq_rep3 | supt16h_3 |


## Fly ([GSE149339](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149336))


  **GEO Status**: Public on May 10, 2020

  **GEO Title**: Pioneer factor GAF cooperates with PBAP and NURF to regulate transcription
  
  **GEO Summary**: The Drosophila pioneer factor GAF is known to be essential for RNA Pol II promoter-proximal pausing and the removal of nucleosomes from a set of target promoters with GAGAG motifs. We and others have speculated that GAF recruits the ISWI family ATP-dependent chromatin remodeling complex NURF, on the basis that NURF and GAF are both required to remodel nucleosomes on an hsp70 promoter in vitro and that GAF interacts physically with NURF. However, GAF was also recently shown to interact with PBAP, a SWI/SNF family remodeler. To test which of these remodeling complexes GAF works with, we depleted GAF, NURF301, BAP170, and NURF301+BAP170 in Drosophila S2 cells using RNAi. We used a combination of PRO-seq, ATAC-seq, 3'RNA-seq, and CUT&RUN to demonstrate that while GAF and PBAP synergistically open chromatin at target promoters which allows Pol II recruitment and pausing to proceed, GAF and NURF also synergistically position the +1 nucleosome to ensure efficient pause release and transition to productive elongation.
  
  **GEO Design**: We treated two independent replicates of Drosophila S2 cells with dsRNA to LACZ (control), GAF, NURF301, BAP170 (the unique subunits of the NURF and PBAP complexes, respectively), and NURF301+BAP170. After 5 days, we harvested cells, validated knockdowns, and performed PRO-seq, ATAC-seq and 3'RNA-seq. We also performed CUT&RUN for both GAF and NURF301 in untreated S2 cells.
  
  **Citation**:  [*Judd, J., Duarte, F. M. & Lis, J. T. Pioneer-like factor GAF cooperates with PBAP (SWI/SNF) and NURF (ISWI) to regulate transcription. Genes Dev. 35, 147–156 (2021).*](http://genesdev.cshlp.org/content/35/1-2/147)

  **Abstract**: Transcriptionally silent genes must be activated throughout development. This requires nucleosomes be removed from promoters and enhancers to allow transcription factor (TF) binding and recruitment of coactivators and RNA polymerase II (Pol II). Specialized pioneer TFs bind nucleosome-wrapped DNA to perform this chromatin opening by mechanisms that remain incompletely understood. Here, we show that GAGA factor (GAF), a Drosophila pioneer-like factor, functions with both SWI/SNF and ISWI family chromatin remodelers to allow recruitment of Pol II and entry to a promoter-proximal paused state, and also to promote Pol II's transition to productive elongation. We found that GAF interacts with PBAP (SWI/SNF) to open chromatin and allow Pol II to be recruited. Importantly, this activity is not dependent on NURF as previously proposed; however, GAF also synergizes with NURF downstream from this process to ensure efficient Pol II pause release and transition to productive elongation, apparently through its role in precisely positioning the +1 nucleosome. These results demonstrate how a single sequence-specific pioneer TF can synergize with remodelers to activate sets of genes. Furthermore, this behavior of remodelers is consistent with findings in yeast and mice, and likely represents general, conserved mechanisms found throughout eukarya. 
  
  **Samples**:
| experiment_accession | run_accession | sample_alias | library_layout | library_strategy |        sample_title        | sample_id  |
|:--------------------:|:-------------:|:------------:|:--------------:|:----------------:|:--------------------------:|:----------:|
|      SRX8174034      |  SRR11607688  |  GSM4498282  |     PAIRED     |     ATAC-seq     |     LACZ_ATACseq_Rep1      |   ctl_1    |
|      SRX8174034      |  SRR11607689  |  GSM4498282  |     PAIRED     |     ATAC-seq     |     LACZ_ATACseq_Rep1      |   ctl_1    |
|      SRX8174035      |  SRR11607690  |  GSM4498283  |     PAIRED     |     ATAC-seq     |     LACZ_ATACseq_Rep2      |   ctl_2    |
|      SRX8174035      |  SRR11607691  |  GSM4498283  |     PAIRED     |     ATAC-seq     |     LACZ_ATACseq_Rep2      |   ctl_2    |
|      SRX8174036      |  SRR11607692  |  GSM4498284  |     PAIRED     |     ATAC-seq     |      GAF_ATACseq_Rep1      |   gaf_1    |
|      SRX8174036      |  SRR11607693  |  GSM4498284  |     PAIRED     |     ATAC-seq     |      GAF_ATACseq_Rep1      |   gaf_1    |
|      SRX8174037      |  SRR11607675  |  GSM4498285  |     PAIRED     |     ATAC-seq     |      GAF_ATACseq_Rep2      |   gaf_2    |
|      SRX8174037      |  SRR11607694  |  GSM4498285  |     PAIRED     |     ATAC-seq     |      GAF_ATACseq_Rep2      |   gaf_2    |
|      SRX8174038      |  SRR11607676  |  GSM4498286  |     PAIRED     |     ATAC-seq     |    BAP170_ATACseq_Rep1     |   b170_1   |
|      SRX8174038      |  SRR11607677  |  GSM4498286  |     PAIRED     |     ATAC-seq     |    BAP170_ATACseq_Rep1     |   b170_1   |
|      SRX8174039      |  SRR11607678  |  GSM4498287  |     PAIRED     |     ATAC-seq     |    BAP170_ATACseq_Rep2     |   b170_2   |
|      SRX8174039      |  SRR11607679  |  GSM4498287  |     PAIRED     |     ATAC-seq     |    BAP170_ATACseq_Rep2     |   b170_2   |
|      SRX8174040      |  SRR11607680  |  GSM4498288  |     PAIRED     |     ATAC-seq     |    NURF301_ATACseq_Rep1    |   n301_1   |
|      SRX8174040      |  SRR11607681  |  GSM4498288  |     PAIRED     |     ATAC-seq     |    NURF301_ATACseq_Rep1    |   n301_1   |
|      SRX8174041      |  SRR11607682  |  GSM4498289  |     PAIRED     |     ATAC-seq     |    NURF301_ATACseq_Rep2    |   n301_2   |
|      SRX8174041      |  SRR11607683  |  GSM4498289  |     PAIRED     |     ATAC-seq     |    NURF301_ATACseq_Rep2    |   n301_2   |
|      SRX8174042      |  SRR11607684  |  GSM4498290  |     PAIRED     |     ATAC-seq     | NURF301BAP170_ATACseq_Rep1 | n301b170_1 |
|      SRX8174042      |  SRR11607685  |  GSM4498290  |     PAIRED     |     ATAC-seq     | NURF301BAP170_ATACseq_Rep1 | n301b170_1 |
|      SRX8174043      |  SRR11607686  |  GSM4498291  |     PAIRED     |     ATAC-seq     | NURF301BAP170_ATACseq_Rep2 | n301b170_2 |
|      SRX8174043      |  SRR11607687  |  GSM4498291  |     PAIRED     |     ATAC-seq     | NURF301BAP170_ATACseq_Rep2 | n301b170_2 |
|      SRX8174044      |  SRR11607698  |  GSM4498295  |     SINGLE     |     RNA-Seq      |      GAF_RNAseq_Rep2       |   gaf_2    |
|      SRX8174045      |  SRR11607699  |  GSM4498296  |     SINGLE     |     RNA-Seq      |     BAP170_RNAseq_Rep1     |   b170_1   |
|      SRX8174046      |  SRR11607700  |  GSM4498297  |     SINGLE     |     RNA-Seq      |     BAP170_RNAseq_Rep2     |   b170_2   |
|      SRX8174047      |  SRR11607701  |  GSM4498298  |     SINGLE     |     RNA-Seq      |    NURF301_RNAseq_Rep1     |   n301_1   |
|      SRX8174048      |  SRR11607702  |  GSM4498299  |     SINGLE     |     RNA-Seq      |    NURF301_RNAseq_Rep2     |   n301_2   |
|      SRX8174049      |  SRR11607703  |  GSM4498300  |     SINGLE     |     RNA-Seq      | NURF301BAP170_RNAseq_Rep1  | n301b170_1 |
|      SRX8174050      |  SRR11607704  |  GSM4498301  |     SINGLE     |     RNA-Seq      | NURF301BAP170_RNAseq_Rep2  | n301b170_2 |
|      SRX8174051      |  SRR11607695  |  GSM4498292  |     SINGLE     |     RNA-Seq      |      LACZ_RNAseq_Rep1      |   ctl_1    |
|      SRX8174052      |  SRR11607696  |  GSM4498293  |     SINGLE     |     RNA-Seq      |      LACZ_RNAseq_Rep2      |   ctl_2    |
|      SRX8174053      |  SRR11607697  |  GSM4498294  |     SINGLE     |     RNA-Seq      |      GAF_RNAseq_Rep1       |   gaf_1    |



## Mouse ([GSE193393](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE193393))


  **GEO Status**: Public on Jun 23, 2022

  **GEO Title**: PHF20 Activates Autophagy Genes through Enhancer Activation via H3K36me2 Binding Activity
  
  **GEO Summary**: Autophagy is a catabolic pathway that maintains cellular homeostasis under various stress conditions, including nutrient-deprived conditions. To elevate autophagic flux to a sufficient level under stress conditions, transcriptional activation of autophagy genes occurs to replenish autophagy components. Here, using combination of RNA-seq, ATAC-seq and ChIP-seq, we demonstrated found that plant homeodomain finger protein 20 (Phf20PHF20), which is an epigenetic reader possessing methyl binding activity, plays a key role in controlling the expression of autophagy genes. PHF20 activates autophagy genes through enhancer activation via H3K36me2 binding activity as an epigenetic reader and that our findings emphasize the importance of nuclear regulation of autophagy.
  
  **GEO Design**: mRNA-seq, ATAC-seq and ChIP-seq experiments under normal and 24hrs of glucose starvation condition in WT and Phf20-/- MEFs
  
  **Citation**: Park SW, Kim J, Oh S, Lee J et al. PHF20 is crucial for epigenetic control of starvation-induced autophagy through enhancer activation. Nucleic Acids Res 2022 Aug 12;50(14):7856-7872. PMID: 35821310
  
  **Samples**:
  
  | experiment_accession | run_accession | sample_alias | library_layout | library_strategy | sample_id  |
  |:--------------------:|:-------------:|:------------:|:--------------:|:----------------:|:----------:|
  |     SRX13654174      |  SRR17483668  |  GSM5776726  |     PAIRED     |     ATAC-seq     |    Wt_1    |
  |     SRX13654175      |  SRR17483667  |  GSM5776727  |     PAIRED     |     ATAC-seq     |    Wt_2    |
  |     SRX13654176      |  SRR17483666  |  GSM5776728  |     PAIRED     |     ATAC-seq     | WtStarv_1  |
  |     SRX13654177      |  SRR17483665  |  GSM5776729  |     PAIRED     |     ATAC-seq     | WtStarv_2  |
  |     SRX13654178      |  SRR17483664  |  GSM5776730  |     PAIRED     |     ATAC-seq     |   Phf_1    |
  |     SRX13654179      |  SRR17483663  |  GSM5776731  |     PAIRED     |     ATAC-seq     |   Phf_2    |
  |     SRX13654180      |  SRR17483662  |  GSM5776732  |     PAIRED     |     ATAC-seq     | PhfStarv_1 |
  |     SRX13654181      |  SRR17483661  |  GSM5776733  |     PAIRED     |     ATAC-seq     | PhfStarv_2 |
  |     SRX13705091      |  SRR17535397  |  GSM5799507  |     PAIRED     |     RNA-Seq      |    Wt_1    |
  |     SRX13705092      |  SRR17535396  |  GSM5799508  |     PAIRED     |     RNA-Seq      |    Wt_2    |
  |     SRX13705093      |  SRR17535395  |  GSM5799509  |     PAIRED     |     RNA-Seq      | WtStarv_1  |
  |     SRX13705094      |  SRR17535394  |  GSM5799510  |     PAIRED     |     RNA-Seq      | WtStarv_2  |
  |     SRX13705095      |  SRR17535392  |  GSM5799511  |     PAIRED     |     RNA-Seq      |   Phf_1    |
  |     SRX13705096      |  SRR17535391  |  GSM5799512  |     PAIRED     |     RNA-Seq      |   Phf_2    |
  |     SRX13705097      |  SRR17535390  |  GSM5799513  |     PAIRED     |     RNA-Seq      | PhfStarv_1 |
  |     SRX13705098      |  SRR17535393  |  GSM5799514  |     PAIRED     |     RNA-Seq      | PhfStarv_2 |
  
  