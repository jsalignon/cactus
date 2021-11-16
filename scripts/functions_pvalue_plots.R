
# ## setting color palettes
# palette_continuous = colorRampPalette(RColorBrewer::brewer.pal(9, 'YlOrRd'), space = 'Lab')(100)
# 
# # discrete palette
# get_pal_disc_name <- function(signed_padj) {
#   if(signed_padj) 'PRGn' else 'Greens'
# }
# pal_disc_dir = -1
# # if(signed_padj) 'PRGn' else 'YlGnBu'
# 
# # continuous palette
# get_pal_cont_name <- function(signed_padj) {
#   if(signed_padj) 'BrBG' else 'Blues'
# }


# col_points_enriched = col_points[6:1]
# col_points_depleted = col_points[11:6]

## getting plot binned and loglog
# add_var_to_plot can take these arguments: none, L2OR, padj_loglog or ov_da

get_plot_binned <- function(p1, signed_padj = F, add_var_to_plot = 'none', add_number = F, point_size = 2, add_number_column = 'ov_da'){
  
  col_bars = rev(RColorBrewer::brewer.pal(11, 'PRGn'))
  col_bars_enriched = col_bars[1:6]

  if(!signed_padj) col_bars = col_bars_enriched
  p_binned = p1 + aes(fill = padj_binned) + scale_fill_manual(values = col_bars, drop = F)
  
  col_points = RColorBrewer::brewer.pal(11, 'BrBG')
  col_points_depleted = col_points[1:6]
  col_points_enriched = col_points[6:11]
  
  add_col = add_var_to_plot
  if(add_col != 'none'){
    if(add_col == 'L2OR' & signed_padj == F) p1$data[[add_col]] %<>% abs
    foo = p1$data[[add_col]]
    all_enriched = all(foo > 0)
    all_depleted = all(foo < 0)
    foo[is.infinite(foo)] = NA
    if(!all(is.na(foo))){
      if(!is.na(all_enriched) & ( all_enriched | !signed_padj ) ) col_points = col_points_enriched
      if(!is.na(all_depleted) & all_depleted) col_points = col_points_depleted
      if(add_col == 'ov_da') col_points = col_points_enriched
      
      if(signed_padj & !all_enriched & !all_depleted & add_col != 'ov_da'){
        # setting white at zero
        rgn = range(foo, na.rm = T)
        colours_break_values = unique(c(seq(rgn[1], 0, len = 6), seq(0, rgn[2], len = 6)))
        my_scale_colour = scale_colour_gradientn(colours = col_points, values = scales::rescale(colours_break_values))
      } else {
        my_scale_colour = scale_colour_gradientn(colours = col_points)
      }
    } else { my_scale_colour = scale_colour_gradientn(colours = col_points) }
    p_binned$data[[add_col]][is.infinite(foo)] = foo
    p_binned = p_binned + geom_point(aes_string(colour = add_col), shape = 19, size = point_size, alpha = 1) + my_scale_colour
  }
  
  if(add_number) {
    nb_cond = length(unique(p1$data$comp_FC))
    rescale_vec = c(nb_cond, seq(0, 20, len = 5))
    overlap_text_size = scales::rescale(rescale_vec, c(6, 3))[1]
    p_binned = p_binned + geom_text(aes_string(label = add_number_column), size = overlap_text_size)
  }
  
  return(p_binned)
}


# pdf('test.pdf')
#   print(p_binned)
# dev.off()

# # # p_binned = get_plot_binned(p1, signed_padj, add_loglog = F, add_L2OR = T)
# p_binned = get_plot_binned(p1, signed_padj, add_loglog, add_L2OR)
# 
# pdf(paste0(key, '.pdf'))
#   print(p_binned)
# dev.off()





get_plot_loglog <- function(p1, add_binned = T, point_size = 0.5, signed_padj = F){
  p_loglog = p1 + aes(fill = padj_loglog) + scale_fill_gradientn(colours = palette_continuous)
  if(add_binned) p_loglog = p_loglog + geom_point(aes(colour = padj_binned), shape = 19, size = point_size, alpha = 1) + scale_color_brewer(palette = get_pal_disc_name(signed_padj), drop = FALSE, direction = pal_disc_dir)
  return(p_loglog)
}


getting_yaxis_terms <- function(df, data_type){
  vec = get_chrom_states_names_vec()
  yaxis_terms = switch(data_type,
    genes = df$term,
    genes_grouped = get_shorter_names(df$term),
    chrom_states = vec[df$chrom_states],
    CHIP = df$CHIP,
    motifs = df$TF
  )
  return(yaxis_terms)
}



getting_padj_loglog_and_binned <- function(df, data_type, signed_padj){
  padj = df$padj
  padj_loglog = get_pval_loglog(padj)
  
  if(signed_padj) {
    signs       = sign(df$L2OR)
    padj        = padj        * signs
    padj_loglog = padj_loglog * signs
  }
  padj_binned = get_padj_binned(padj, data_type, signed_padj)
  
  df$padj_loglog = padj_loglog
  df$padj_binned = padj_binned
  
  return(df)
}


## getting pval columns
get_padj_binned_breaks <- function(data_type){
  breaks = switch(data_type,
    genes_self     = c(1,    0.2, 0.05, 1e-5 , 1e-20 , 1e-100, 0),
    peaks_self     = c(1,    0.2, 0.05, 1e-5 , 1e-20 , 1e-100, 0),
    func_anno      = c(1,    0.2, 0.05, 1e-5 , 1e-20 , 1e-100, 0),
    chrom_states   = c(1,    0.2, 0.05, 1e-5 , 1e-20 , 1e-100, 0),
    CHIP           = c(1,    0.2, 0.05, 1e-5 , 1e-20 , 1e-100, 0),
    motifs         = c(1,    0.2, 0.05, 1e-5 , 1e-20 , 1e-100, 0)
  )
  return(breaks)
}


get_padj_binned <- function(padj, data_type, signed_padj){
  breaks = get_padj_binned_breaks(data_type)
  br_len = length(breaks) - 1
  if(signed_padj) breaks = unique(c(breaks, -breaks))
  padj_binned = cut(padj, breaks = breaks, include.lowest = T, right = F)
  if(signed_padj) {
    lpb = levels(padj_binned)
    not_signif = '> 0.05'
    padj_binned %<>% factor(., levels = c(lpb, not_signif))
    padj_binned[padj_binned %in% lpb[c(1, length(lpb))]] = not_signif
    new_levels_order = c(br_len:2, br_len * 2 + 1, (2 * br_len - 1):(br_len + 1))
    new_levels = rev(levels(padj_binned)[new_levels_order])
    padj_binned = factor(padj_binned, levels = new_levels)
  } 
  return(padj_binned)
}



# get_signed_pval <- function(pval, effect){
#   sign = ifelse(effect > 1, 1, -1)
#   new_pval = pval * sign
#   return(new_pval)
# }

# these functions return negative log

get_pval_log <- function(pval){
  pval_log = -log10(pval)
  pval_log[pval == 0] = 330
  return(pval_log)
}

get_raw_pval_from_pval_log <- function(pval_log){
  pval = 10^-pval_log
  pval[pval_log == 330] = 0
  return(pval)
}

get_pval_loglog <- function(pval){
  pval_log = get_pval_log(pval)
  pval_loglog = log10(pval_log + 1)
  return(pval_loglog)
}

get_raw_pval_from_pval_loglog <- function(pval_loglog){
  pval_log = 10^pval_loglog - 1
  pval = get_raw_pval_from_pval_log(pval_log)
  return(pval)
}

## examples
# vec = c(1, 0, 0.05, 0.001, 10^-10, 10^-300, 10^-400)
# df = cbind(pv = vec, pv_log = get_pval_log(vec), pv_loglog = get_pval_loglog(vec), pv_from_log = get_raw_pval_from_pval_log(get_pval_log(vec)), pv_from_loglog = get_raw_pval_from_pval_loglog(get_pval_loglog(vec)))
# df
#          pv    pv_log pv_loglog pv_from_log pv_from_loglog
# [1,]  1e+00   0.00000 0.0000000       1e+00          1e+00
# [2,]  0e+00 330.00000 2.5198280       0e+00          0e+00
# [3,]  5e-02   1.30103 0.3619223       5e-02          5e-02
# [4,]  1e-03   3.00000 0.6020600       1e-03          1e-03
# [5,]  1e-10  10.00000 1.0413927       1e-10          1e-10
# [6,] 1e-300 300.00000 2.4785665      1e-300         1e-300
# [7,]  0e+00 330.00000 2.5198280       0e+00          0e+00


 ## Notes

 # // note: we convert to log or loglog scale for clustering otherwise small values would be ignored as compared to large one during the clustering step (i.e. pvalues 1 and 9 would have much more weight than 10^-50 and 10^-100 which would both be considered as zero. On a loglog scale the opposite happens. We want a good separation on the low pvalues!).

 # // note: for the loglog we need to ad +1 in the second log call, otherwise pvalues that are 1 would become infinite
 #
 # // # we add 330 since it is the limit of precision of the machine
 # // # the limit of precision on this machine is 10^-323
 # // 10^-324 # 0
 # // 10^-323 # 9.881313e-324
 # // 10^-323 == 0 # F
 # // 10^-324 == 0 # T
 # // -log10(10^-323) # 323
 #


 ## solving the equation (not valid now since use minus log instead of log)
 #           pval_loglog = -log10(-pval_log + 1)
 #          -pval_loglog = log10(-pval_log + 1)
 #       10^-pval_loglog = -pval_log + 1
 #  10^-pval_loglog - 1  = -pval_log
 # -10^-pval_loglog + 1  =  pval_log
