
# this function merges different columns in a single column by collapsing them with a "|"
get_merged_columns <- function(res, sel_columns_chr, char_to_clean = 'PF', collapse = '|'){
  df = res[, sel_columns_chr, drop = F]
  merged_columns = apply(df, 1, which)
  merged_columns_clean = purrr::map_chr(merged_columns, ~paste0(names(df)[.x], collapse = collapse) %>% gsub(paste0(char_to_clean, '_'), '', .)) %>% unname
  return(merged_columns_clean)
}
# FDR_split = c(1.3, 2, 5)
# res = res_simple
# sel_columns_chr =  paste0('FDR_', FDR_split)
# char_to_clean = 'FDR'

filter_entries_by_threshold <- function(res, threshold_type, threshold_value){
  if (threshold_type == 'FDR'){
    res_simple$padj < 10^-threshold_value
  } else if (threshold_type == 'rank'){
    res_simple$rank <= threshold_value
  } else {
    stop('wrong threshold type')
  }
}
