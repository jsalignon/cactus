

get_formatted_table <- function(df){
  nms = names(df)
  if('L2OR' %in% nms) df %<>% dplyr::arrange(padj, desc(abs(L2OR)))
  if('L2FC' %in% nms) df %<>% dplyr::arrange(padj, desc(abs(L2FC)))
  
  # reformatting table
  v_to_round = c('L2FC', 'L2OR')
  for(col in v_to_round){
    if(col %in% nms) df[[col]] %<>% round(2)
  }
  v_to_format = c('pval', 'padj')
  for(col in v_to_format){
    if(col %in% nms) df[[col]] %<>% format(digits = 2)
  }
  
  return(df)
}

