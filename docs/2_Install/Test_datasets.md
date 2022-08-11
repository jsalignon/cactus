

* [Introduction](/README.md): [Quick Start](/docs/1_Intro/Quick_start.md), [Flowchart](/docs/1_Intro/Flowchart.md), [Outputs structure](/docs/1_Intro/Outputs_structure.md)
* [Install](/docs/2_Install.md): [Dependencies](/docs/2_Install/Dependencies.md), [Containers](/docs/2_Install/Containers.md), [Data](/docs/2_Install/Data.md), [Test_datasets](/docs/2_Install/Test_datasets.md)
* [Run](/docs/3_Run/3_Run.md): [Input Data](/docs/3_Run/Input_data.md), [Input Files](/docs/3_Run/Input_files.md), [Parameters](/docs/3_Run/Parameters.md)


[](END_OF_MENU)





## Introduction

Tests datasets have beed generated for each of the 4 species supported by Cactus. This allows to make sure the pipeline works equally well for each species. Downloading the smallest dataset (fly) should be sufficient for users to get a sense on how the pipeline works. Here are the sizes of the test datasets:

| specie | raw_files | sampled_files | sampled_atac | sampled_mrna |
|:------:|:---------:|:-------------:|:------------:|:------------:|
|  fly   |   17 GB   |    377 MB     |    329 MB    |    48 MB     |
|  worm  |   68 GB   |    1.4 GB     |    1.3 GB    |    41 MB     |
| human  |  118 GB   |    3.7 GB     |    3.6 GB    |    190 MB    |
| mouse  |  104 GB   |     12 GB     |    11 GB     |    395 MB    |

And below are details on the test datasets origins and labels.

... add a table and a not on genome and transcriptome size here to explain the differences!



## Worm and human ([GSE98758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98758))

  **GEO Status**: Public on Aug 29, 2018

  **GEO Title**: Genome-wide DNA accessibility maps and differential gene expression using ChIP-seq, ATAC-seq and RNA-seq for the human secondary fibroblast cell line hiF-T and whole worms with and without knockdown of FACT complex
  
  **GEO Summary**: To assess the mechanisms by which FACT depletion leads to increased sensitivity of cells to be reprogrammed, we measured the chromatin accessibility landscape using ATAC-seq following mock treatment, SSRP1 knockdown, or SUPT16H knockdown in human fibroblasts and mock, hmg-3 or hmg-4 knockdown in whole worms, and differential gene expression in hmg-3 knockout mutants or following mock, hmg-4, or spt-16 knockdown by RNAseq.
  
  **GEO Design**: Examination of two FACT complex components in human cells and worms with ChIP-seq, ATAC-seq and RNA-seq
  
  **Citation**: [*Kolundzic E, Ofenbauer A, Bulut SI, Uyar B et al. FACT Sets a Barrier for Cell Fate Reprogramming in Caenorhabditis elegans and Human Cells. Dev Cell 2018 Sep 10;46(5):611-626.e12. PMID: 30078731*](https://www.sciencedirect.com/science/article/pii/S1534580718305598)

  **Abstract**: The chromatin regulator FACT (facilitates chromatin transcription) is essential for ensuring stable gene expression by promoting transcription. In a genetic screen using Caenorhabditis elegans, we identified that FACT maintains cell identities and acts as a barrier for transcription factor-mediated cell fate reprogramming. Strikingly, FACT’s role as a barrier to cell fate conversion is conserved in humans as we show that FACT depletion enhances reprogramming of fibroblasts. Such activity is unexpected because FACT is known as a positive regulator of gene expression, and previously described reprogramming barriers typically repress gene expression. While FACT depletion in human fibroblasts results in decreased expression of many genes, a number of FACT-occupied genes, including reprogramming-promoting factors, show increased expression upon FACT depletion, suggesting a repressive function of FACT. Our findings identify FACT as a cellular reprogramming barrier in C. elegans and humans, revealing an evolutionarily conserved mechanism for cell fate protection.
  
  **Homologs: worms -> humans**
    - hmg-4, hmg-3  -> SSRP1
    -        spt-16 -> SUPT16H

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



## Mouse ([GSE181797](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE181797))


  **GEO Status**: Public on Oct 01, 2021

  **GEO Title**: Revealing regulatory transcription factors of aging in the kidney and kidney through multi-omics analysis
  
  **GEO Summary**: Aging is accompanied by a gradual loss of systemic metabolic homeostasis, which is maintained by multiple-organs, especially the kidney and kidney. However, a systematic study of the regulatory networks and regulatory transcription factors (TFs) of aging in the kidney and kidney remains lacking. Here, we performed an integrated analysis of multi-omics datasets in the kidney and kidney from young and aged mice, including RNA-seq, reduced representation bisulfite sequencing (RRBS) and ATAC-seq datasets, which indicated that enhanced inflammation and dysregulated metabolism were conserved signatures in aged kidney and kidney in both the transcriptome and epigenome. Transcription factor and regulatory network analysis indicated that activation of AP-1 and SPI1 was responsible for enhanced inflammation, and down-regulation of HNFs and PPARs contributed to dysregulated metabolism in aged kidney and kidney. Importantly, we found that the activation of AP-1 was cellular autonomous in aged hepatic and renal cells. However, enhanced SPI1 was caused by elevated infiltration of macrophages. Importantly, inhibition of AP-1 with small molecule combination attenuated inflammation phenotypes of aging in vivo and in vitro. Taken together, our analysis revealed common signatures and regulatory TFs of aging in the kidney and kidney, providing novel targets for the development of anti-aging interventions.
  
  **GEO Design**: mRNA and ATAC-seq profiles of young and old kidney and kidney.
  
  **Citation**:  *Citation missing*
  
  **Samples**:
| experiment_accession | run_accession | sample_alias | library_layout | library_strategy |      sample_title      | sample_id |
|:--------------------:|:-------------:|:------------:|:--------------:|:----------------:|:----------------------:|:---------:|
|     SRX11708663      |  SRR15406500  |  GSM5511490  |     PAIRED     |     RNA-Seq      | Young_kidney_mRNA_rep1 | YngKid_1  |
|     SRX11708664      |  SRR15406501  |  GSM5511491  |     PAIRED     |     RNA-Seq      | Young_kidney_mRNA_rep2 | YngKid_2  |
|     SRX11708665      |  SRR15406502  |  GSM5511492  |     PAIRED     |     RNA-Seq      | Young_kidney_mRNA_rep3 | YngKid_3  |
|     SRX11708666      |  SRR15406503  |  GSM5511493  |     PAIRED     |     RNA-Seq      | Young_kidney_mRNA_rep4 | YngKid_4  |
|     SRX11708667      |  SRR15406504  |  GSM5511494  |     PAIRED     |     RNA-Seq      |  Old_kidney_mRNA_rep1  | OldKid_1  |
|     SRX11708668      |  SRR15406505  |  GSM5511495  |     PAIRED     |     RNA-Seq      |  Old_kidney_mRNA_rep2  | OldKid_2  |
|     SRX11708669      |  SRR15406506  |  GSM5511496  |     PAIRED     |     RNA-Seq      |  Old_kidney_mRNA_rep3  | OldKid_3  |
|     SRX11708670      |  SRR15406507  |  GSM5511497  |     PAIRED     |     RNA-Seq      |  Old_kidney_mRNA_rep4  | OldKid_4  |
|     SRX11708671      |  SRR15406508  |  GSM5511498  |     PAIRED     |     RNA-Seq      | Young_liver_mRNA_rep1  | YngLiv_1  |
|     SRX11708672      |  SRR15406509  |  GSM5511499  |     PAIRED     |     RNA-Seq      | Young_liver_mRNA_rep2  | YngLiv_2  |
|     SRX11708673      |  SRR15406510  |  GSM5511500  |     PAIRED     |     RNA-Seq      | Young_liver_mRNA_rep3  | YngLiv_3  |
|     SRX11708674      |  SRR15406511  |  GSM5511501  |     PAIRED     |     RNA-Seq      | Young_liver_mRNA_rep4  | YngLiv_4  |
|     SRX11708675      |  SRR15406512  |  GSM5511502  |     PAIRED     |     RNA-Seq      |  Old_liver_mRNA_rep1   | OldLiv_1  |
|     SRX11708676      |  SRR15406513  |  GSM5511503  |     PAIRED     |     RNA-Seq      |  Old_liver_mRNA_rep2   | OldLiv_2  |
|     SRX11708677      |  SRR15406514  |  GSM5511504  |     PAIRED     |     RNA-Seq      |  Old_liver_mRNA_rep3   | OldLiv_3  |
|     SRX11708678      |  SRR15406515  |  GSM5511505  |     PAIRED     |     RNA-Seq      |  Old_liver_mRNA_rep4   | OldLiv_4  |
|     SRX11708679      |  SRR15406516  |  GSM5511506  |     PAIRED     |     ATAC-seq     | Young_kidney_ATAC_rep1 | YngKid_1  |
|     SRX11708680      |  SRR15406517  |  GSM5511507  |     PAIRED     |     ATAC-seq     | Young_kidney_ATAC_rep2 | YngKid_2  |
|     SRX11708681      |  SRR15406518  |  GSM5511508  |     PAIRED     |     ATAC-seq     | Young_kidney_ATAC_rep3 | YngKid_3  |
|     SRX11708682      |  SRR15406519  |  GSM5511509  |     PAIRED     |     ATAC-seq     |  Old_kidney_ATAC_rep1  | OldKid_1  |
|     SRX11708683      |  SRR15406520  |  GSM5511510  |     PAIRED     |     ATAC-seq     |  Old_kidney_ATAC_rep2  | OldKid_2  |
|     SRX11708684      |  SRR15406521  |  GSM5511511  |     PAIRED     |     ATAC-seq     |  Old_kidney_ATAC_rep3  | OldKid_3  |
|     SRX11708685      |  SRR15406522  |  GSM5511512  |     PAIRED     |     ATAC-seq     | Young_liver_ATAC_rep1  | YngLiv_1  |
|     SRX11708686      |  SRR15406523  |  GSM5511513  |     PAIRED     |     ATAC-seq     | Young_liver_ATAC_rep2  | YngLiv_2  |
|     SRX11708687      |  SRR15406524  |  GSM5511514  |     PAIRED     |     ATAC-seq     | Young_liver_ATAC_rep3  | YngLiv_3  |
|     SRX11708688      |  SRR15406525  |  GSM5511515  |     PAIRED     |     ATAC-seq     |  Old_liver_ATAC_rep1   | OldLiv_1  |
|     SRX11708689      |  SRR15406526  |  GSM5511516  |     PAIRED     |     ATAC-seq     |  Old_liver_ATAC_rep2   | OldLiv_2  |
|     SRX11708690      |  SRR15406527  |  GSM5511517  |     PAIRED     |     ATAC-seq     |  Old_liver_ATAC_rep3   | OldLiv_3  |
