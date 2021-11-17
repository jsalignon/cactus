
################################################
## Ordering the x-axis terms (comparisons)

# example of input: comp_order = 'g4d7_vs_g4d3|g4d3_vs_g4d1|g4d7_vs_g4d1' ; up_down_pattern = 'UDUD'
# up_down_pattern: can be either UDUD (up, down, up, down...) or UUDD (up, up, ..., down, down ...) 

get_comp_order_levels <- function(comp_order, up_down_pattern = 'UUDD'){
  comp_order1 = gsub('_vs_', '_', comp_order)
  comp_order1 = strsplit(comp_order1, '|', fixed = T)[[1]]
  len = length(comp_order1)
  if(up_down_pattern == 'UUDD') { 
    times_UD  = 1 ; each_UD   = len
    times_CFC = 2   ; each__CFC = 1
  }
  if(up_down_pattern == 'UDUD') { 
    times_UD  = len ; each_UD   = 1
    times_CFC = 1 ; each__CFC = 2
  }
  up_down     = rep(c('up','down'), times = times_UD,  each = each_UD)
  comp_order2 = rep(comp_order1   , times = times_CFC, each = each__CFC)
  comp_order3 = paste0(comp_order2, '_', up_down)
  
  return(comp_order3)
}



################################################
## Selecting the y-axis terms (enriched gene sets, motifs...)

# the input should be a matrix with rownames containing the terms of interest

# the output is the selected terms ordered by clustering

select_y_axis_terms_grouped_plot <- function(mat, nshared = 6, nunique = 20, ntotal = 26, threshold_type = 'fixed', threshold_value = 0.05, seed = 38, remove_similar = F, remove_similar_n = 2, reverse = T){

  set.seed(seed)

  if(reverse) mat = -mat
  mat_full = mat

  # if there are less terms than we want to select we downscale
  if( nrow(mat) < ntotal ){

    mat_final = mat_full

  } else {

    if(remove_similar){

      # keeping only the top n yaxis_terms for yaxis_terms with similar names
      rn = rownames(mat)
      yaxis_terms_names = purrr::map_chr(strsplit(rn, '_'), 1)
      duplicated_yaxis_terms = names(which(table(yaxis_terms_names) > 1))
      worst_dupli = c()
      for(yaxis_terms in duplicated_yaxis_terms){
        sel = which(yaxis_terms_names == yaxis_terms)
        dupli_min_pval = apply(mat[sel, ], 1, min)
        ranked = rank(dupli_min_pval, ties.method = 'random')
        worst_dupli1 = names(ranked)[ranked > remove_similar_n]
        worst_dupli = c(worst_dupli, worst_dupli1)
      }
      mat %<>% .[!rownames(.) %in% worst_dupli, ]

    }

    ## selecting the top 6 shared terms
    all_pvalues = c(unname(mat))
    if(threshold_type == 'quantile') threshold = quantile(all_pvalues, threshold_value)
    if(threshold_type == 'fixed') threshold = threshold_value
    shared_sum = apply(mat, 1, function(x) sum(x < threshold))
    shared_terms = c()
    if(any(shared_sum > 1)){
      ranked = rank(shared_sum, ties.method = 'random')
      shared_terms = names(ranked)[ranked <= nshared]
      # we remove selected shared terms so that they do not appear in the unique terms
      mat %<>% .[!rownames(.) %in% shared_terms, ]
    }

    # selecting the top_N terms for each comparison
    ncomp = ncol(mat)
    top_N = floor(nunique / ncomp)
    top_N_terms = apply(mat, 2, function(x) { ranked = rank(x, ties.method = 'random') ; names(ranked)[ranked <= top_N] } ) %>% c %>% unique
    mat %<>% .[!rownames(.) %in% top_N_terms, ]

    # filling with the lowest pvalues overall
    cur_terms = unique(c(shared_terms, top_N_terms))
    nmissing = ntotal - length(cur_terms)
    if(nmissing > 0){
      ranked = rank(apply(mat, 1, min), ties.method = 'random')
      low_pval_terms = names(ranked)[ranked <= nmissing]
      cur_terms %<>% c(., low_pval_terms)
    }

    mat_final = mat_full %>% .[rownames(.) %in% cur_terms, ]
  }

  # finally clustering terms
  col_order = hclust(dist(as.matrix(mat_final)))$order
  terms_levels = rownames(mat_final)[col_order]

  return(terms_levels)

}


## debugging:

# gene enrichment terms
# nshared = 6 ; nunique = 20 ; ntotal = 26 ; threshold_type = 'fixed'; threshold_value = 0.05 ; seed = 38

# Transcription Factors
# nshared = 8 ; nunique = 25 ; ntotal = 40 ; threshold_type = 'quantile'; threshold_value = 0.25 ; seed = 38





################################################
## Clustering, reformating to long format, and adding padj columns

# The input of this function should be a matrix of padj_loglog values that is already ordered (with a pre-set order, or via clustering with the  select_y_axis_terms_grouped_plot function).
# The output is dataframe in long format to serve as input to ggplot2

add_matrix_indexes_to_df <- function(cur_mat, cur_df, rows, cols, data_type, signed_loglog = F){

  # conversion to long format
  df1 = data.frame(yaxis_terms = rownames(cur_mat), cur_mat, stringsAsFactors = F)
  names(df1)[-1] =  paste0('padj_loglog.', names(df1)[-1])
  df2 = reshape(df1, direction = 'long', timevar = 'comp_FC', varying = names(df1)[-1])
  rownames(df2) <- df2$id <- df2$padj_loglog <- NULL

  # adding heatmap row and column indexes
  df2$rowInd = rep(1:rows, times = cols)
  df2$colInd = rep(1:cols, each = rows)
  
  # adding these new columns to the df
  df3 = dplyr::inner_join(cur_df, df2, by = c('yaxis_terms', 'comp_FC'))

  return(df3)
}


################################################
## Plotting the heatmap

getting_heatmap_base <- function(cur_df, rows, cols, title, cur_mat){
  p1 = ggplot(data = cur_df, aes(x = colInd - 0.5, y = rowInd - 0.5, xmin = colInd - 1, xmax = colInd, ymin = rowInd - 1, ymax = rowInd)) + geom_rect() + theme_bw() + scale_x_continuous(breaks = (1:cols) - 0.5, labels = colnames(cur_mat), expand = c(0, 0)) + scale_y_continuous(breaks = (1:rows) - 0.5, labels = rownames(cur_mat), expand = c(0, 0)) + coord_fixed(1) + theme(axis.title = element_blank(), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 8)) + ggtitle(title)
  return(p1)
}



# getting_heatmap_base <- function(cur_df, rows, cols, title, cur_mat){
#   p1 = ggplot(data = cur_df, 
#     aes(x = colInd - 0.5, y = rowInd - 0.5, xmin = colInd - 1, xmax = colInd, ymin = rowInd - 1, ymax = rowInd)) + 
#     geom_rect() + theme_bw() + 
#     scale_x_continuous(breaks = (1:cols) - 0.5, labels = colnames(cur_mat), expand = c(0, 0), na.value="white" ) + 
#     scale_y_continuous(breaks = (1:rows) - 0.5, labels = rownames(cur_mat), expand = c(0, 0)) + 
#     coord_fixed(1) 
#     + theme(axis.title = element_blank(), axis.text.y = element_text(size = 7), axis.text.x = element_text(size = 7, angle = 90, hjust = 1, vjust = 0.5), axis.line = element_line(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), plot.title = element_text(hjust = 0.5, size = 8)) + 
#     ggtitle(title)
#   return(p1)
# }
# 
# 













