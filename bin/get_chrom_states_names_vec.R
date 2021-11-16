
## terms for chromatin states analysis. See reference below
## hiHMM: Bayesian non-parametric joint inference of chromatin state maps
# https://www.nature.com/articles/nature13415
# https://www.nature.com/articles/nature13415/figures/2

get_chrom_states_names_vec <- function(){
  vec = c('1. Promoter', '2. Enhancer 1', '3. Enhancer 2', '4. Transcription 5\' 1', '5. Transcription 5\' 2', '6. Gene, H420me1', '7. Transcription 3\' 1', '8. Transcription 3\' 2', '9. Transcription 3\' 3', '10. PC repressed 1', '11. PC repressed 2', '12. Heterochromatin 1', '13. Heterochromatin 2', '14. Low signal 1', '15. Low signal 2', '16. Low signal 3') %>% setNames(., c('1_Pro', '2_Enh1', '3_Enh2', '4_Egn1', '5_Egn2', '6_Egn3', '7_Egn4', '8_Egn5', '9_Egn6', '10_Rep1', '11_Rep2', '12_Het1', '13_Het2', '14_Low1', '15_Low2', '16_Low3'))
  return(vec)
}
