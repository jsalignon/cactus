
get_rank_per_unique_character <- function(v_char){
  for(char in unique(v_char)){
    sel = which(v_char == char)
    v_char[sel] = seq_len(length(sel))
  }
  v_char
}


get_new_name_by_unique_character <- function(v_char, clean = T){
	u_char = get_rank_per_unique_character(v_char)
	v_char = paste(v_char, u_char, sep = '_')
  if(clean) v_char = gsub('_1$', '', v_char)
	v_char
}

get_shorter_names <- function(long_name, max_char = 40){
  short_name = substr(long_name, 1, max_char)
  short_name %<>% gsub(' $', '', .)
  
  # if shortening has merged go terms, then we add an unique id
  has_merged = duplicated(long_name) != duplicated(short_name)
  if(any(has_merged)){
    has_merged1 = short_name %in% short_name[has_merged]
    short_name[has_merged1] %<>% get_new_name_by_unique_character
  } 
  
  return(short_name)
}





# get_comp_order_levels <- function(comp_order){
#   comp_order1 = gsub('_vs_', '_', comp_order)
#   comp_order1 = strsplit(comp_order1, '|', fixed = T)[[1]]
#   up_down = rep(c('up','down'), length(comp_order1))
#   comp_order1 = rep(comp_order1, each = 2)
#   comp_order1 = paste0(comp_order1, '_', up_down)
#   return(comp_order1)
# }

# example of input: comp_order = 'g4d7_vs_g4d3|g4d3_vs_g4d1|g4d7_vs_g4d1'


# my_palette <- function(x) colorRampPalette(rev(RColorBrewer::brewer.pal(9, 'YlOrRd')), space='Lab')(x)
# 
# 
# add_pval_loglog_and_binned <- function(df){
#   df$padj_loglog = -log(-log(df$padj))
#   breaks = c(0.05, 1e-2, 1e-3, 1e-5, 1e-10, 1e-25, 1e-50, 0)
#   df$padj_binned = cut(df$padj, breaks = breaks)
#   return(df)
# }
# 
# # input of this function should be a df from a GRanges object with these columns in order: "chr, start, end, name, score, strand", and optional fields after
# 
# export_df_to_bed <- function(df, filename){
#   df$start = df$start - 1
#   df$start[df$start < 0] = 0
#   write.table(df, file = filename, quote = F, sep = '\t', row.names = F, col.names = F)
# }
# 
# 
# get_prom_bed_df_table <- function(promoters_df, mRNA_res_table){
#   prom_bed_df = promoters_df
#   prom_bed_df %<>% dplyr::inner_join(mRNA_res_table[, c('gene_id', 'pval')], by = 'gene_id')
#   prom_bed_df %<>% dplyr::mutate(score = round(-log10(pval), 2), name = gene_name)
#   prom_bed_df %<>% dplyr::select(chr, start, end, name, score, strand, gene_id)
#   return(prom_bed_df)
# }



