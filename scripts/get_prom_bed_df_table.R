

get_prom_bed_df_table <- function(promoters_df, mRNA_res_table){
  prom_bed_df = promoters_df
  prom_bed_df %<>% dplyr::inner_join(mRNA_res_table[, c('gene_id', 'pval')], by = 'gene_id')
  prom_bed_df %<>% dplyr::mutate(score = round(-log10(pval), 2), name = gene_name)
  prom_bed_df %<>% dplyr::select(chr, start, end, name, score, strand, gene_id)
  return(prom_bed_df)
}



