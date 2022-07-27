
## Nature 17 states HiHMM model 
# states details 1: https://www.nature.com/articles/nature13415/figures/2
# states details 2: https://gander.wustl.edu/cgi-bin/hgTables?db=dm6&hgta_track=hiHMM_models&hgta_table=hiHMM_M1K16_L3&hgta_doSchema=describe+table+schema
#       reference: https://www.nature.com/articles/nature13415

## ENCODE 18 states ChromHMM model 
#              states details: https://hzhou.scholar.harvard.edu/blog/chromhmm
#  66 mouse samples reference: https://www.nature.com/articles/s42003-021-01756-4
# 833 human samples reference: https://www.nature.com/articles/s41586-020-03145-z

## checking that we have indeed the right number of ENCODE files
# cd ${cactus}/data
# ls -p mouse/chromatin_states | grep -v / | wc -l # 66
# ls -p human/chromatin_states | grep -v / | wc -l # 833

get_chrom_states_names_vec <- function(){
  vec = c(
    
    # HiHMM
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
    '17_Unmap' = '17. Unmappable',
    
    # ENCODE
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
    
  return(vec)
}

# R
# library(magrittr)
# v_hihmm = list.files('human/chromatin_states/iHMM.M1K16.human_H1') %>% gsub('.bed', '', .) 
# v_chromhmm = list.files('human/chromatin_states/ENCFF306VHM') %>% gsub('.bed', '', .) 
# get_chrom_states_names_vec()[c(v_hihmm, v_chromhmm)]

# v_hihmm = list.files('human/chromatin_states/iHMM.M1K16.human_H1') %>% gsub('.bed', '', .) %>% {gtools::mixedsort(.)} %>% writeLines
# list.files('human/chromatin_states/ENCFF306VHM') %>% gsub('.bed', '', .) %>% {gtools::mixedsort(.)} %>% writeLines





# get_chrom_states_names_vec <- function(){
#   vec = c('1. Promoter', '2. Enhancer 1', '3. Enhancer 2', '4. Transcription 5\' 1', '5. Transcription 5\' 2', '6. Gene, H420me1', '7. Transcription 3\' 1', '8. Transcription 3\' 2', '9. Transcription 3\' 3', '10. PC repressed 1', '11. PC repressed 2', '12. Heterochromatin 1', '13. Heterochromatin 2', '14. Low signal 1', '15. Low signal 2', '16. Low signal 3') %>% setNames(., c('1_Pro', '2_Enh1', '3_Enh2', '4_Egn1', '5_Egn2', '6_Egn3', '7_Egn4', '8_Egn5', '9_Egn6', '10_Rep1', '11_Rep2', '12_Het1', '13_Het2', '14_Low1', '15_Low2', '16_Low3'))
#   return(vec)
# }

