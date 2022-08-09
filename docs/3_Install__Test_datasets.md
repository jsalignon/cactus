
* Introduction: [Cactus](/README.md), [Quick Start](/docs/1_Intro__Quick_start.md), 
* Overview: [Graph](/docs/2_Overview__Graph.md), [Outputs structure](/docs/2_Overview__Outputs_structure.md)
* Install: [Dependencies](/docs/3_Install__Dependencies.md), [Containers](/docs/3_Install__Containers.md), [Data](/docs/3_Install__Data.md), [Test_datasets](/docs/3_Install__Test_datasets.md)
* Run: [Input Data](/docs/4_Run__Input_data.md), [Input Files](/docs/4_Run__Input_files.md), [Parameters](/docs/4_Run__Parameters.md)
* Preprocessing: ATAC: [Method](/docs/5_AP__Method.md), [Figures](/docs/5_AP__Figures.md), [MultiQC](/docs/5_AP__MultiQC.md), mRNA: [Method](/docs/6_MP__Method.md), [MultiQC](/docs/6_MP__MultiQC.md)
* Differential Abundance: [ATAC](/docs/7_DA__DiffBind.md), [mRNA](/docs/7_DA__Sleuth.md), [Figures](/docs/7_DA__Figures.md), [Tables](/docs/7_DA__Tables.md)
* Splitting peak sets: [Split](/docs/8_SP__Splitting_peak_sets.md), [Venn Diagrams](/docs/8_SP__Venn_diagrams.md)
* Enrichment: [ChIP](/docs/9_Enrich__CHIP.md), [Chromatin state](/docs/9_Enrich__Chromatin_states.md), [Motifs](/docs/9_Enrich__Motifs.md), [Func. Anno.](/docs/9_Enrich__Functional_annotations.md), [Figures](/docs/9_Enrich__Figures.md), [Tables](/docs/9_Enrich__Tables.md)


[](END_OF_MENU)


## Worm and human ([GSE98758](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE98758))

  **GEO Status**: Public on Aug 29, 2018

  **GEO Title**: Genome-wide DNA accessibility maps and differential gene expression using ChIP-seq, ATAC-seq and RNA-seq for the human secondary fibroblast cell line hiF-T and whole worms with and without knockdown of FACT complex
  
  **GEO Summary**: To assess the mechanisms by which FACT depletion leads to increased sensitivity of cells to be reprogrammed, we measured the chromatin accessibility landscape using ATAC-seq following mock treatment, SSRP1 knockdown, or SUPT16H knockdown in human fibroblasts and mock, hmg-3 or hmg-4 knockdown in whole worms, and differential gene expression in hmg-3 knockout mutants or following mock, hmg-4, or spt-16 knockdown by RNAseq.
  
  **GEO Design**: Examination of two FACT complex components in human cells and worms with ChIP-seq, ATAC-seq and RNA-seq
  
  **Citation**: [*Kolundzic E, Ofenbauer A, Bulut SI, Uyar B et al. FACT Sets a Barrier for Cell Fate Reprogramming in Caenorhabditis elegans and Human Cells. Dev Cell 2018 Sep 10;46(5):611-626.e12. PMID: 30078731*](https://www.sciencedirect.com/science/article/pii/S1534580718305598)

  **abstract**: The chromatin regulator FACT (facilitates chromatin transcription) is essential for ensuring stable gene expression by promoting transcription. In a genetic screen using Caenorhabditis elegans, we identified that FACT maintains cell identities and acts as a barrier for transcription factor-mediated cell fate reprogramming. Strikingly, FACT’s role as a barrier to cell fate conversion is conserved in humans as we show that FACT depletion enhances reprogramming of fibroblasts. Such activity is unexpected because FACT is known as a positive regulator of gene expression, and previously described reprogramming barriers typically repress gene expression. While FACT depletion in human fibroblasts results in decreased expression of many genes, a number of FACT-occupied genes, including reprogramming-promoting factors, show increased expression upon FACT depletion, suggesting a repressive function of FACT. Our findings identify FACT as a cellular reprogramming barrier in C. elegans and humans, revealing an evolutionarily conserved mechanism for cell fate protection.
  
  **homologs: worms -> humans**
    - hmg-4, hmg-3  -> SSRP1
    -        spt-16 -> SUPT16H


## Fly ([GSE149339](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149336))


  **GEO Status**: Public on May 10, 2020

  **GEO Title**: Pioneer factor GAF cooperates with PBAP and NURF to regulate transcription
  
  **GEO Summary**: The Drosophila pioneer factor GAF is known to be essential for RNA Pol II promoter-proximal pausing and the removal of nucleosomes from a set of target promoters with GAGAG motifs. We and others have speculated that GAF recruits the ISWI family ATP-dependent chromatin remodeling complex NURF, on the basis that NURF and GAF are both required to remodel nucleosomes on an hsp70 promoter in vitro and that GAF interacts physically with NURF. However, GAF was also recently shown to interact with PBAP, a SWI/SNF family remodeler. To test which of these remodeling complexes GAF works with, we depleted GAF, NURF301, BAP170, and NURF301+BAP170 in Drosophila S2 cells using RNAi. We used a combination of PRO-seq, ATAC-seq, 3'RNA-seq, and CUT&RUN to demonstrate that while GAF and PBAP synergistically open chromatin at target promoters which allows Pol II recruitment and pausing to proceed, GAF and NURF also synergistically position the +1 nucleosome to ensure efficient pause release and transition to productive elongation.
  
  **GEO Design**: We treated two independent replicates of Drosophila S2 cells with dsRNA to LACZ (control), GAF, NURF301, BAP170 (the unique subunits of the NURF and PBAP complexes, respectively), and NURF301+BAP170. After 5 days, we harvested cells, validated knockdowns, and performed PRO-seq, ATAC-seq and 3'RNA-seq. We also performed CUT&RUN for both GAF and NURF301 in untreated S2 cells.
  
  **Citation**:  [*Judd, J., Duarte, F. M. & Lis, J. T. Pioneer-like factor GAF cooperates with PBAP (SWI/SNF) and NURF (ISWI) to regulate transcription. Genes Dev. 35, 147–156 (2021).*](http://genesdev.cshlp.org/content/35/1-2/147)

  **abstract**: Transcriptionally silent genes must be activated throughout development. This requires nucleosomes be removed from promoters and enhancers to allow transcription factor (TF) binding and recruitment of coactivators and RNA polymerase II (Pol II). Specialized pioneer TFs bind nucleosome-wrapped DNA to perform this chromatin opening by mechanisms that remain incompletely understood. Here, we show that GAGA factor (GAF), a Drosophila pioneer-like factor, functions with both SWI/SNF and ISWI family chromatin remodelers to allow recruitment of Pol II and entry to a promoter-proximal paused state, and also to promote Pol II's transition to productive elongation. We found that GAF interacts with PBAP (SWI/SNF) to open chromatin and allow Pol II to be recruited. Importantly, this activity is not dependent on NURF as previously proposed; however, GAF also synergizes with NURF downstream from this process to ensure efficient Pol II pause release and transition to productive elongation, apparently through its role in precisely positioning the +1 nucleosome. These results demonstrate how a single sequence-specific pioneer TF can synergize with remodelers to activate sets of genes. Furthermore, this behavior of remodelers is consistent with findings in yeast and mice, and likely represents general, conserved mechanisms found throughout eukarya. 
  

## Mouse ([GSE149339](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149336))

=> use this paper instead!
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE162048


    ############################################################
    ## Fly (GSE149339)
    ########################
    
    Status 	Public on May 10, 2020
    Title 	Pioneer factor GAF cooperates with PBAP and NURF to regulate transcription [ATAC-seq]
    Organism 	Drosophila melanogaster
    Experiment type 	Genome binding/occupancy profiling by high throughput sequencing
    Summary 	The Drosophila pioneer factor GAF is known to be essential for RNA Pol II promoter-proximal pausing and the removal of nucleosomes from a set of target promoters with GAGAG motifs. We and others have speculated that GAF recruits the ISWI family ATP-dependent chromatin remodeling complex NURF, on the basis that NURF and GAF are both required to remodel nucleosomes on an hsp70 promoter in vitro and that GAF interacts physically with NURF. However, GAF was also recently shown to interact with PBAP, a SWI/SNF family remodeler. To test which of these remodeling complexes GAF works with, we depleted GAF, NURF301, BAP170, and NURF301+BAP170 in Drosophila S2 cells using RNAi. We used a combination of PRO-seq, ATAC-seq, 3'RNA-seq, and CUT&RUN to demonstrate that while GAF and PBAP synergistically open chromatin at target promoters which allows Pol II recruitment and pausing to proceed, GAF and NURF also synergistically position the +1 nucleosome to ensure efficient pause release and transition to productive elongation.
      
    Overall design 	We treated two independent replicates of Drosophila S2 cells with dsRNA to LACZ (control), GAF, NURF301, BAP170 (the unique subunits of the NURF and PBAP complexes, respectively), and NURF301+BAP170. After 5 days, we harvested cells, validated knockdowns, and performed PRO-seq, ATAC-seq and 3'RNA-seq. We also performed CUT&RUN for both GAF and NURF301 in untreated S2 cells.
    Web link 	https://doi.org/10.1101/2020.05.10.087262
      
    Contributor(s) 	Judd J, Duarte FM, Lis JT
    Citation missing 	Has this study been published? Please login to update or notify GEO.


      

############################################################
## Mice (GSE181797)
###################

# => no publication yet!

Title 	Revealing regulatory transcription factors of aging in the kidney and kidney through multi-omics analysis

Summary 	Aging is accompanied by a gradual loss of systemic metabolic homeostasis, which is maintained by multiple-organs, especially the kidney and kidney. However, a systematic study of the regulatory networks and regulatory transcription factors (TFs) of aging in the kidney and kidney remains lacking. Here, we performed an integrated analysis of multi-omics datasets in the kidney and kidney from young and aged mice, including RNA-seq, reduced representation bisulfite sequencing (RRBS) and ATAC-seq datasets, which indicated that enhanced inflammation and dysregulated metabolism were conserved signatures in aged kidney and kidney in both the transcriptome and epigenome. Transcription factor and regulatory network analysis indicated that activation of AP-1 and SPI1 was responsible for enhanced inflammation, and down-regulation of HNFs and PPARs contributed to dysregulated metabolism in aged kidney and kidney. Importantly, we found that the activation of AP-1 was cellular autonomous in aged hepatic and renal cells. However, enhanced SPI1 was caused by elevated infiltration of macrophages. Importantly, inhibition of AP-1 with small molecule combination attenuated inflammation phenotypes of aging in vivo and in vitro. Taken together, our analysis revealed common signatures and regulatory TFs of aging in the kidney and kidney, providing novel targets for the development of anti-aging interventions.




############################################################
## Fly (GSE149339)
########################

Status 	Public on May 10, 2020
Title 	Pioneer factor GAF cooperates with PBAP and NURF to regulate transcription [ATAC-seq]
Organism 	Drosophila melanogaster
Experiment type 	Genome binding/occupancy profiling by high throughput sequencing
Summary 	The Drosophila pioneer factor GAF is known to be essential for RNA Pol II promoter-proximal pausing and the removal of nucleosomes from a set of target promoters with GAGAG motifs. We and others have speculated that GAF recruits the ISWI family ATP-dependent chromatin remodeling complex NURF, on the basis that NURF and GAF are both required to remodel nucleosomes on an hsp70 promoter in vitro and that GAF interacts physically with NURF. However, GAF was also recently shown to interact with PBAP, a SWI/SNF family remodeler. To test which of these remodeling complexes GAF works with, we depleted GAF, NURF301, BAP170, and NURF301+BAP170 in Drosophila S2 cells using RNAi. We used a combination of PRO-seq, ATAC-seq, 3'RNA-seq, and CUT&RUN to demonstrate that while GAF and PBAP synergistically open chromatin at target promoters which allows Pol II recruitment and pausing to proceed, GAF and NURF also synergistically position the +1 nucleosome to ensure efficient pause release and transition to productive elongation.
  
Overall design 	We treated two independent replicates of Drosophila S2 cells with dsRNA to LACZ (control), GAF, NURF301, BAP170 (the unique subunits of the NURF and PBAP complexes, respectively), and NURF301+BAP170. After 5 days, we harvested cells, validated knockdowns, and performed PRO-seq, ATAC-seq and 3'RNA-seq. We also performed CUT&RUN for both GAF and NURF301 in untreated S2 cells.
Web link 	https://doi.org/10.1101/2020.05.10.087262
  
Contributor(s) 	Judd J, Duarte FM, Lis JT
Citation missing 	Has this study been published? Please login to update or notify GEO.
