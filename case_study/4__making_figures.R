
# homedir=~
# eval homedir=$homedir
# cactus_dir=$homedir/workspace/cactus
# cd $cactus_dir/case_study/figures


R


# loading functions and libraries

library(ggplot2)
library(magrittr)
library(data.table)

readRDS1 <- function(x) readRDS(paste0('panels/', x))
pdf_a4   <- function(x, ...) pdf(x, paper = 'a4r', ...)
update_text_size <- function(px) {
	px + theme(
		axis.text.y = element_text(size = 8), 
		axis.text.x = element_text(size = 10), 
		plot.title = element_text(size = 10, hjust = 1)
		)
}

# Fig 5: peak_self heatmaps

p1 = readRDS1('mRNA__Null__1.3__ctl__peaks_self__heatmap__worm.rds') + ggtitle('worm, DEGs') 
p2 = readRDS1('mRNA__Null__1.3__ctl__peaks_self__heatmap__human.rds') + ggtitle('human, DEGs') 
p3 = readRDS1('ATAC__all__1.3__ctl__peaks_self__heatmap__worm.rds') + ggtitle('worm, DARs') 
p4 = readRDS1('ATAC__all__1.3__ctl__peaks_self__heatmap__human.rds') + ggtitle('human, DARs') 

pp = ggpubr::ggarrange(p1, p2, p3, p4, ncol = 2, nrow = 2, common.legend = T, 
	legend = 'right', labels = letters[1:4])

pdf('Figure_5.pdf', width = 7.5, height = 7.5)
	print(pp)
dev.off()


# Fig 6: CHIP and motifs heatmaps

p1 = readRDS1('ATAC__all__1.3__ctl__motifs__heatmap__worm.rds') + ggtitle('worm, motifs') 
p2 = readRDS1('ATAC__all__1.3__ctl__CHIP__heatmap__worm.rds') + ggtitle('worm, CHIP') 
p3 = readRDS1('ATAC__all__1.3__ctl__motifs__heatmap__human.rds') + ggtitle('human, motifs') 
p4 = readRDS1('ATAC__all__1.3__ctl__CHIP__heatmap__human.rds') + ggtitle('human, CHIP') 

p1 %<>% update_text_size 
p2 %<>% update_text_size 
p3 %<>% update_text_size 
p4 %<>% update_text_size 

p1 = p1 + theme(
	legend.text = element_text(size = 6), 
	legend.title = element_text(size = 6),
	legend.key.size = unit(0.5, 'cm'))

pp = ggpubr::ggarrange(p1, p2, p3, p4, ncol = 4, nrow = 1, common.legend = T, 
	legend = 'right', labels = letters[1:4], widths = c(1.19, 1.3, 1.02, 1.4))

pdf('Figure_6.pdf', width = 7, height = 5.7)
	print(pp)
dev.off()



# Fig 7: CHIP both_atac and both_mrna heatmaps

p1 = readRDS1('both_ATAC__all__1.3__ctl__CHIP__heatmap.rds') + ggtitle('human, DARs, CHIP,\n co-regulated genes')
p2 = readRDS1('both_mRNA__all__1.3__ctl__CHIP__heatmap.rds') + ggtitle('human, promoters, CHIP\n co-regulated genes')
p1 = p1 + theme(plot.title = element_text(hjust = 1))
p2 = p2 + theme(plot.title = element_text(hjust = 1))

pp = ggpubr::ggarrange(p1, p2, ncol = 2, nrow = 1, common.legend = T, 
	legend = 'right', labels = letters[1:2])

pdf('Figure_7.pdf', width = 7, height = 7)
	print(pp)
dev.off()


# Fig 8: chrom states heatmaps

p1 = readRDS1('ATAC__all__1.3__ctl__chrom_states__heatmap__worm.rds') + ggtitle('HiHMM, worm\nL3 larvae')
p2 = readRDS1('ATAC__all__1.3__ctl__chrom_states__heatmap__human_2.rds') + ggtitle('ChromHMM, human\nfibroblast cell line')
p3 = readRDS1('ATAC__all__1.3__ctl__chrom_states__heatmap__human_1.rds') + ggtitle('HiHMM, human\nlymphoblastoid cell line')

p1 %<>% update_text_size 
p2 %<>% update_text_size 
p3 %<>% update_text_size 

p1 = p1 + theme(
	legend.text = element_text(size = 6), 
	legend.title = element_text(size = 6),
	legend.key.size = unit(0.5, 'cm')
	)

pp = ggpubr::ggarrange(p1, p2, p3, ncol = 3, nrow = 1, common.legend = T, 
	legend = 'right', labels = letters[1:3], widths = c(1.2, 1.03, 1.27))

pdf('Figure_8.pdf', width = 7, height = 3.5)
	print(pp)
dev.off()



# Tab S1: reprogramming genes defined in Figure 6L of the Kolundzic et al. study

dt = openxlsx::read.xlsx('panels/res_simple.xlsx') %>% setDT
dt = dt[ET == 'mRNA' & COMP %in% c('ssrp1_vs_ctl', 'supt16h_vs_ctl')]

# reprogramming-promoting genes
dt_up = dt[gene_name %in% c('FGF2', 'CEBPB', 'ESRRB', 'PRDM15', 'BMP2', 'JMJD2C', 'SALL4', 'LIN28B', 'KDM4C')]
dt_up = dt_up[order(-L2FC), c('COMP', 'gene_name', 'gene_id', 'pval', 'padj', 'L2FC', 'FC')]

# reprogramming-preventing genes
dt_dw = dt[gene_name %in% c('SUV39H1', 'SUV39H2', 'PTPN11', 'NR2F1', 'SUMO2', 'CHAF1B', 'RBBP7', 'DRAM1')]
dt_dw = dt_dw[order(-L2FC), c('COMP', 'gene_name', 'gene_id', 'pval', 'padj', 'L2FC', 'FC')]

openxlsx::write.xlsx(list(up = dt_up, down = dt_dw), 'Table_S1.xlsx')



homedir=~
eval homedir=$homedir
cactus_dir=$homedir/workspace/cactus
cd $cactus_dir/case_study/figures

# ~/workspace/cactus/docs/examples/png
# /raid/jersal/workspace/cactus/docs/examples/png

library(ggplot2)

png_to_gg <- function(png_file, margin = 20) {
	ggplot2::ggplot() + ggplot2::annotation_custom(
		grid::rasterGrob(png::readPNG(png_file), 
			width = ggplot2::unit(1,"npc"), 
			height = ggplot2::unit(1,"npc")), 
		-Inf, Inf, -Inf, Inf) + theme(plot.margin = 
    margin(t = margin, r = margin, b = margin, l = margin, unit = "pt"))
}

plot_by_chunk <- function(lp, nrow, ncol){
	nplots = length(lp)
	npages = ceiling(nplots / (nrow * ncol))
	lp1 = vector('list', length = npages)

	for(page in 1:npages){
		cur_panel = 1 + (page - 1) * 12
		v_panels = cur_panel:min(cur_panel + 11, nplots)
		p1 = ggpubr::ggarrange(plotlist = lp[v_panels], nrow = 4, ncol = 3, 
			labels = v_panels)
		lp1[[page]] = p1
	}

	return(lp1)
}


v_files = list.files('../../docs/examples/png', full.names = T)

v_files_1 = c(
	grep('\\/ctl_1__', v_files, value = T),
	grep('without_control_pca', v_files, value = T),
	grep('spearman_correlation_heatmap', v_files, value = T),
	grep('ATAC__peaks__grouped', v_files, value = T),
	grep('hmg4_vs_ctl__ATAC_', v_files, value = T),
	grep('hmg4_vs_ctl__mRNA_', v_files, value = T),
	# grep('mRNA_volcano__', v_files, value = T),
	grep('venn_up', v_files, value = T),
	grep('__barplot', v_files, value = T),
	grep('__heatmap\\.', v_files, value = T)
	)
v_files %>% .[!. %in% v_files_1]
# [1] "../../docs/examples/png/hmg4_vs_spt16__ATAC_volcano__no_rtr.png"
# [2] "../../docs/examples/png/hmg4_vs_spt16__ATAC_volcano.png"
# [3] "../../docs/examples/png/mRNA_volcano__no_gtr.png"
# [4] "../../docs/examples/png/mRNA_volcano__with_gtr.png"

v_files_1 %>% .[which(duplicated(.))] # character(0)

v_files_1 %<>% .[. != '../../docs/examples/png/hmg4_vs_ctl__ATAC_other_plots.png']
v_files_1 %<>% .[. != '../../docs/examples/png/hmg4_vs_ctl__mRNA_other_plots.png']

lp = lapply(v_files_1, png_to_gg)

lp1 = plot_by_chunk(lp, 4, 3)

pdf('Sup_fig_test.pdf', width = 7, height = 10)
	print(lp1)
dev.off()





# lp = lapply(v_files_1, png_to_gg)
# p1 = ggpubr::ggarrange(plotlist = lp, nrow = 4, ncol = 3, 
# 	labels = 'auto')

# pdf('Sup_fig_test.pdf', width = 7, height = 10)
# 	print(p1)
# dev.off()

