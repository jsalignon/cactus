
## terms for chromatin states analysis. See reference below
## hiHMM: Bayesian non-parametric joint inference of chromatin state maps
# https://www.nature.com/articles/nature13415
# https://www.nature.com/articles/nature13415/figures/2

## details on the 18 states model from ENCODE ChromHMM chromatin state files
# https://hzhou.scholar.harvard.edu/blog/chromhmm


get_chrom_states_names_vec <- function(){
  vec = c(
    
    # HiHMM
    '1_Pro'   = '1. Promoter',
    '2_Enh1'  = '2. Enhancer 1',
    '3_Enh2'  = '3. Enhancer 2',
    '4_Egn1'  = '4. Transcription 5\' 1',
    '5_Egn2'  = '5. Transcription 5\' 2',
    '6_Egn3'  = '6. Gene, H420me1',
    '7_Egn4'  = '7. Transcription 3\' 1',
    '8_Egn5'  = '8. Transcription 3\' 2',
    '9_Egn6'  = '9. Transcription 3\' 3',
    '10_Rep1' = '10. PC repressed 1',
    '11_Rep2' = '11. PC repressed 2',
    '12_Het1' = '12. Heterochromatin 1',
    '13_Het2' = '13. Heterochromatin 2',
    '14_Low1' = '14. Low signal 1',
    '15_Low2' = '15. Low signal 2',
    '16_Low3' = '16. Low signal 3',
    
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

