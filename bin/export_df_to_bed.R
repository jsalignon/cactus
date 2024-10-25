
# input of this function should be a df from a GRanges object with these columns in order: "chr, start, end, name, score, strand", and optional fields after

export_df_to_bed <- function(df, filename){
  df$chr = paste0('chr', df$chr)
  df$start = df$start - 1
  df$start[df$start < 0] = 0
  df$start %<>% as.integer
  df$end   %<>% as.integer
  write.table(df, file = filename, quote = F, sep = '\t', row.names = F, col.names = F)
}
