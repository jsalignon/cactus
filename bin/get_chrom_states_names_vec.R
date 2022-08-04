
## Nature 17 states HiHMM model 
# states details 1: https://www.nature.com/articles/nature13415/figures/2
# states details 2: https://gander.wustl.edu/cgi-bin/hgTables?db=dm6&hgta_track=hiHMM_models&hgta_table=hiHMM_M1K16_L3&hgta_doSchema=describe+table+schema
#       reference: https://www.nature.com/articles/nature13415

## ENCODE 18 states ChromHMM model 
#        states details human: https://hzhou.scholar.harvard.edu/blog/chromhmm
#        states details mouse: https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs42003-021-01756-4/MediaObjects/42003_2021_1756_Fig1_HTML.png?as=webp
#  66 mouse samples reference: https://www.nature.com/articles/s42003-021-01756-4
# 833 human samples reference: https://www.nature.com/articles/s41586-020-03145-z

## checking that we have indeed the right number of ENCODE files
# cd ${cactus}/data
# ls -p mouse/chromatin_states | grep -v / | wc -l # 66
# ls -p human/chromatin_states | grep -v / | wc -l # 833


# HiHMM
chrom_state_worm <- function(){
  c(
    '1_Pro'    = '1. Promoter',
    '2_Enh1'   = '2. Enhancer 1',
    '3_Enh2'   = '3. Enhancer 2',
    '4_Egn1'   = '4. Transcription 5\' 1',
    '5_Egn2'   = '5. Transcription 5\' 2',
    '6_Egn3'   = '6. Gene, H420me1',
    '7_Egn4'   = '7. Transcription 3\' 1',
    '8_Egn5'   = '8. Transcription 3\' 2',
    '9_Egn6'   = '9. Transcription 3\' 3',
    '10_Rep1'  = '10. PC repressed 1',
    '11_Rep2'  = '11. PC repressed 2',
    '12_Het1'  = '12. Heterochromatin 1',
    '13_Het2'  = '13. Heterochromatin 2',
    '14_Low1'  = '14. Low signal 1',
    '15_Low2'  = '15. Low signal 2',
    '16_Low3'  = '16. Low signal 3',
    '17_Unmap' = '17. Unmappable'
  )
}

chrom_state_human <- function(){
  c(
    TssA     = '1. Active TSS',
    TssFlnk  = '2. Flanking TSS',
    TssFlnkU = '3. Flanking TSS Upstream',
    TssFlnkD = '4. Flanking TSS Downstream',
    Tx       = '5. Strong transcription',
    TxWk     = '6. Weak transcription',
    EnhG1    = '7. Genic enhancer1',
    EnhG2    = '8. Genic enhancer2',
    EnhA1    = '9. Active Enhancer 1',
    EnhA2    = '10. Active Enhancer 2',
    EnhWk    = '11. Weak Enhancer',
    ZNF_Rpts = '12. ZNF genes & repeats',
    Het      = '13. Heterochromatin',
    TssBiv   = '14. Bivalent/Poised TSS',
    EnhBiv   = '15. Bivalent Enhancer',
    ReprPC   = '16. Repressed PolyComb',
    ReprPCWk = '17. Weak Repressed PolyComb',
    Quies    = '18. Quiescent/Low'
  )
}

chrom_state_mouse <- function(){
  c(
    Tss      = '1. Active TSS',
    TssFlnk  = '2. Flanking TSS',
    Tx       = '3. Transcription',
    TxWk     = '4. Weak transcription',
    EnhG     = '5. Enhancer in gene',
    Enh      = '6. Enhancer',
    EnhLo    = '7. Weak enhancer',
    EnhPois  = '8. Poised enhancer',
    EnhPr    = '9. Primed enhancer',
    TssBiv   = '10. Bivalent TSS',
    ReprPC   = '11. Repressed by PolyComb',
    ReprPCWk = '12. Repressed by PolyComb (wk)',
    QuiesG   = '13. Quiescent gene',
    Quies    = '14. Quiescent',
    Quies2   = '15. Quiescent',
    Quies3   = '16. Quiescent',
    Quies4   = '17. Quiescent',
    Het      = '18. Heterochromatin'
  )
}



get_chrom_states_names_vec <- function(raw_names){
  
  if('1_Pro'    %in% raw_names) vec = chrom_state_worm()
  if('TssFlnkU' %in% raw_names) vec = chrom_state_human()
  if('Quies4'   %in% raw_names) vec = chrom_state_mouse()
    
  return(vec)
  
}

# R
# library(magrittr)
# v_hihmm = list.files('human/chromatin_states/iHMM.M1K16.human_H1') %>% gsub('.bed', '', .) 
# v_chromhmm = list.files('human/chromatin_states/ENCFF306VHM') %>% gsub('.bed', '', .) 
# get_chrom_states_names_vec()[c(v_hihmm, v_chromhmm)]

# v_hihmm = list.files('human/chromatin_states/iHMM.M1K16.human_H1') %>% gsub('.bed', '', .) %>% {gtools::mixedsort(.)} %>% writeLines
# list.files('human/chromatin_states/ENCFF306VHM') %>% gsub('.bed', '', .) %>% {gtools::mixedsort(.)} %>% writeLines


# ls data/mouse/chromatin_states/*/
# ls data/human/chromatin_states/*/

# data/mouse/chromatin_states/ENCFF940GBJ/:
# Enh.bed   EnhLo.bed    EnhPr.bed  Quies2.bed  Quies4.bed  QuiesG.bed  ReprPCWk.bed  TssBiv.bed   Tx.bed
# EnhG.bed  EnhPois.bed  Het.bed    Quies3.bed  Quies.bed   ReprPC.bed  Tss.bed       TssFlnk.bed  TxWk.bed

# data/human/chromatin_states/ENCFF770KKC/:
# EnhA1.bed  EnhBiv.bed  EnhG2.bed  Het.bed    ReprPC.bed    TssA.bed    TssFlnk.bed   TssFlnkU.bed  TxWk.bed
# EnhA2.bed  EnhG1.bed   EnhWk.bed  Quies.bed  ReprPCWk.bed  TssBiv.bed  TssFlnkD.bed  Tx.bed        ZNF_Rpts.bed

