

plot_venn_diagrams <- function(lgenes, prefix){
  
  stsp = strsplit(names(lgenes), '_')
  FCs = sapply(stsp, function(x) x[1])
  ETs = sapply(stsp, function(x) x[2])
  
  line_colours = c(up = '#800080', down = '#008000')
  area_colours = c(atac = '#0b5d95', mrna = '#ff6600')
  
  lg = length(lgenes)
  two_sets = lg == 2
  
  
  col = line_colours[ETs]
  fill = area_colours[FCs]
  if(two_sets){
    len = sapply(lgenes, length)
    rotation_degree = ifelse(len[1] < len[2], 180, 0)
    euler_d = T 
    scaled = T
    tl = 3 # potential text items length
    filename = paste0(prefix, '__venn_up_or_down.pdf')
  } else {
    rotation_degree = 0
    euler_d = F 
    scaled = F
    tl = 15
    filename = paste0(prefix, '__venn_up_and_down.pdf')
  }
  
  v1 <- venn.diagram(lgenes, filename = NULL, lwd = 4, fill = fill, col = col, ext.text = FALSE, rotation.degree = rotation_degree, euler.d = euler_d, scaled = scaled, alpha = rep(0.5, lg), fontfamily = rep('sans', tl), cex = rep(1.6, tl))
  # , label.col = rep("white", lg), 
  gh <- grobHeight(v1)
  
  if(two_sets){
    # vjust = ( max(len) / min(len) / 5 ) - 1
    # vjust = ifelse(vjust > 12, 12, vjust)
    # vjust = ifelse(vjust < 0, 0, vjust)
    vjust = 12 # + 4 * ( log10(min(len)) -   log10(max(len)) )
  } else {
    vjust = 0
  }
  
  y0 = unit(0.5, 'npc') + 0.7 * gh[1]
  title <- textGrob(prefix, y = y0, vjust = vjust, gp = gpar(fontsize = 10, fontfamily = 'sans'))
  
  pdf(filename)
    grid.newpage()
    grid.draw(v1)
    grid.draw(title)
  dev.off()

}

# y0 = unit(0.5, 'npc') 
# title <- textGrob(prefix, y = y0, gp = gpar(fontsize = 10, fontfamily = 'sans'))


## this would be the solution to solve the title problem. However, gridExtra and ggplot2 are not available in this container...
# pdf(filename)
#   grid::grid.arrange(v1, v1, ncol=2, main = "Main title")
#   # grid.newpage()
#   # grid.draw(v1)
#   # grid.draw(title)
# dev.off()

## some testing code
# PF1 = 'all'
# FDR1 = 1.3
# atac_up = get_file_name('ATAC', PF1, 'up', FDR1)
# atac_down = get_file_name('ATAC', PF1, 'down', FDR1)
# mrna_up = get_file_name('mRNA', 'Null', 'up', FDR1)
# mrna_down = get_file_name('mRNA', 'Null', 'down', FDR1)
# lgenes_two_sets_up = list(atac_up = readRDS(atac_up), mrna_up = readRDS(mrna_up))
# lgenes_two_sets_down = list(atac_down = readRDS(atac_down), mrna_down = readRDS(mrna_down))
# lgenes_four_sets = list(atac_up = readRDS(atac_up), mrna_up = readRDS(mrna_up), atac_down = readRDS(atac_down), mrna_down = readRDS(mrna_down))
# plot_venn_diagrams(lgenes_two_sets_down, 'test_two_down')
# plot_venn_diagrams(lgenes_two_sets_up, 'test_two_up')
# plot_venn_diagrams(lgenes_four_sets, 'test_four')

# lgenes_two_sets_up1 = list(atac_up = c(readRDS(mrna_up)[1:100], readRDS(atac_up)), mrna_up = readRDS(mrna_up))
# plot_venn_diagrams(lgenes_two_sets_up1, 'test_two_up1')

